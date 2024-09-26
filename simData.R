library(phytools)
library(geiger)
library(tidyverse)
library(cowplot)
library(mvMORPH)
library(RevGadgets)

fitContinuousMCMC <- function(
    phy,
    d,
    model = c("BM", "Trend", "SSP", "ACDC.exp", "ACDC.lin", "lambda"),
    Ngens = 1000000,
    sampleFreq = 1000,
    printFreq = 1000,
    propwidth = rep(0.1, 5),
    node.priors = NULL,
    root.prior = NULL,
    acdc.prior = NULL,
    sample.node.states = TRUE,
    outputName = "mcmc.output"
) {
  
  # Supporting Functions
  .vcv.anc <- function(tree, anc.nodes = TRUE) {
    n <- length(tree$tip.label)
    h <- .nodeheights(tree)[order(tree$edge[, 2]), 2]
    h <- c(h[1:n], 0, h[(n + 1):length(h)])
    M <- mrca(tree, full = anc.nodes)[c(1:n, anc.nodes * (n + 2:tree$Nnode)), c(1:n, anc.nodes * (n + 2:tree$Nnode))]
    C <- matrix(h[M], nrow(M), ncol(M))
    if (anc.nodes)
      rownames(C) <- colnames(C) <- c(tree$tip.label, n + 2:tree$Nnode)
    else rownames(C) <- colnames(C) <- tree$tip.label
    return(C)
  }
  
  .nodeheights <- function(tree) {
    if (attr(tree, "order") != "cladewise" || is.null(attr(tree, "order")))
      t <- reorder(tree)
    else t <- tree
    root <- length(t$tip) + 1
    X <- matrix(NA, nrow(t$edge), 2)
    for (i in 1:nrow(t$edge)) {
      if (t$edge[i, 1] == root) {
        X[i, 1] <- 0
        X[i, 2] <- t$edge.length[i]
      } else {
        X[i, 1] <- X[match(t$edge[i, 1], t$edge[, 2]), 2]
        X[i, 2] <- X[i, 1] + t$edge.length[i]
      }
    }
    if (attr(tree, "order") != "cladewise" || is.null(attr(tree, "order")))
      o <- apply(matrix(tree$edge[, 2]), 1, function(x, y) which(x == y), y = t$edge[, 2])
    else o <- 1:nrow(t$edge)
    return(X[o, ])
  }
  
  transform_lambda <- function(phy, lambda) {
    phy$edge.length <- phy$edge.length * lambda
    return(phy)
  }
  
  .slater.getvcv <- function(C, model, params) {
    if (model == "BM" | model == "Trend") {
      V <- C * params$sigma
    }
    
    if (model == "ACDC.exp") {
      if (!params$scaler == 0) {
        C <- (exp(params$scaler * C) - 1) / params$scaler
      }
      V <- C * params$sigma
    }
    
    if (model == "ACDC.lin") {
      V <-  params$sigma * (C + ((params$scaler * (C^2)) / 2))
    }
    
    if (model == "OU" | model == "SSP") {
      oux = .ou.vcv(C)
      if (params$scaler == 0) {
        V <- C * params$sigma
      } else {
        V <- oux(params$scaler)
        V <- params$sigma * V
      }
    }
    
    if (model == "lambda") {
      lambda_phy <- transform_lambda(phy, params$scaler)
      V <- vcv(lambda_phy) * params$sigma
    }
    
    return(V)
  }
  
  .slater.getproposal <- function(currentState, psi, bound = FALSE, min = NA, max = NA, prop.type = "sw") {
    if(prop.type == "sw") {
      prop <- currentState + ((runif(1) - 0.5) * psi)
    } else if(prop.type == "mlt") {
      u <- runif(1)
      prop <- currentState * exp((u-0.5)*psi)
    }
    
    if(bound == TRUE) {
      if(is.na(min) & is.na(max)) stop("you must specify limits if using a bounded proposal!")
      if(!is.na(min) & is.na(max)) max <- Inf
      if(is.na(min) & !is.na(max)) min <- -Inf
      if (prop < min) { prop <- min + (min - prop) } # if lower than lower bound, bounce back
      if (prop > max) { prop <- max - (prop - max) } # if higher than upper bound, bounce back
    }
    
    if(prop.type == "sw") return(prop)
    if(prop.type == "mlt") return(list(prop = prop, hastings = u))
  }
  
  .slater.generateproposal <- function(ngen, propwidth, model, params, node.priors, root.prior, scaler.prior = NULL, sampling.nodes) {
    pwRoot <- propwidth[1]
    pwRate <- propwidth[2]
    pwScale <- propwidth[3]
    pwTheta <- propwidth[4]
    pwFossil <- propwidth[5]
    
    if(!is.null(node.priors) | sampling.nodes == TRUE) {
      if(model == "BM") {
        if(ngen %% 3 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 3 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE, min = -Inf, max = Inf))
        }
        if(ngen %% 3 == 0) {
          fu <- .slater.updatewhichnode(model, ngen, params$node.states)
          params$node.states[fu] <- .slater.nodeupdate(fu, model, params$node.states, node.priors, pwFossil)
        }
      } else if (model == "Trend") {
        if(ngen %% 4 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 4 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
        }
        if(ngen %% 4 == 3) {
          params$scaler <- .slater.getproposal(params$scaler, pwScale, bound = FALSE)
        }
        if(ngen %% 4 == 0) {
          fu <- .slater.updatewhichnode(model, ngen, params$node.states)
          params$node.states[fu] <- .slater.nodeupdate(fu, model, params$node.states, node.priors, pwFossil)
        }
      } else if(model =="ACDC.lin" | model == "ACDC.exp") {
        if(ngen %% 3 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 3 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
          params$scaler <- .slater.getproposal(params$scaler, pwScale, bound = TRUE, min = scaler.prior$min, max = scaler.prior$max)
        }
        if(ngen %% 3 == 0) {
          fu <- .slater.updatewhichnode(model, ngen, params$node.states)
          params$node.states[fu] <- .slater.nodeupdate(fu, model, params$node.states, node.priors, pwFossil)
        }
      } else if (model == "SSP") {
        if(ngen %% 4 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 4 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
        }
        if(ngen %% 4 == 3) {
          params$scaler <- exp(.slater.getproposal(log(params$scaler), pwScale, bound = TRUE, min = log(10^-5), max = Inf))
        }
        if(ngen %% 4 == 0) {
          fu <- .slater.updatewhichnode(model, ngen, params$node.states)
          params$node.states[fu] <- .slater.nodeupdate(fu, model, params$node.states, node.priors, pwFossil)
        }
      } else if(model == "OU") {
        if(ngen %% 5 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 5 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
        }
        if(ngen %% 5 == 3) {
          params$scaler <- exp(.slater.getproposal(log(params$scaler), pwScale, bound = TRUE, min =  log(10^-5), max = Inf))
        }
        if(ngen %% 5 == 4) {
          params$theta <- .slater.getproposal(params$theta, pwTheta)
        }
        if(ngen %% 5 == 0) {
          fu <- .slater.updatewhichnode(model, ngen, params$node.states)
          params$node.states[fu] <- .slater.nodeupdate(fu, model, params$node.states, node.priors, pwFossil)
        }
      }
    } else  {
      if(model == "BM") {
        if(ngen %% 2 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 2 == 0) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE, min = -Inf, max = Inf))
        }
      } else if (model == "Trend") {
        if(ngen %% 3 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 3 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
        }
        if(ngen %% 3 == 0) {
          params$scaler <- .slater.getproposal(params$scaler, pwScale, bound = FALSE)
        }
      } else if(model =="ACDC.lin" | model == "ACDC.exp") {
        if(ngen %% 2 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 2 == 0) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
          params$scaler <- .slater.getproposal(params$scaler, pwScale, bound = TRUE, min = scaler.prior$min, max = scaler.prior$max)
        }
      } else if (model == "SSP") {
        if(ngen %% 3 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 3 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
        }
        if(ngen %% 3 == 0) {
          params$scaler <- exp(.slater.getproposal(log(params$scaler), pwScale, bound = TRUE, min = log(10^-5), max = Inf))
        }
      } else if(model == "OU") {
        if(ngen %% 4 == 1) {
          if(is.null(root.prior)) {
            params$root <- .slater.getproposal(params$root, pwRoot)
          } else {
            params$root <- .slater.nodeupdate(1, model, params$root, root.prior, pwRoot)
          }
        }
        if(ngen %% 4 == 2) {
          params$sigma <- exp(.slater.getproposal(log(params$sigma), pwRate, bound= TRUE,  min = -Inf, max = Inf))
        }
        if(ngen %% 4 == 3) {
          params$scaler <- exp(.slater.getproposal(log(params$scaler), pwScale, bound = TRUE, min =  log(10^-5), max = Inf))
        }
        if(ngen %% 4 == 0) {
          params$theta <- .slater.getproposal(params$theta, pwTheta)
        }
      }
    }
    return(params)
  }
  
  .slater.generatestartingvalues <- function(phy, d, C, model, node.priors, root.prior, scaler.prior, sampling.nodes) {
    node.pics <- ace(d, phy, type = "continuous", method = "pic")$ace
    
    if(sampling.nodes == TRUE) {
      node.states <- node.pics[-1]
      if(!is.null(node.priors)) {
        repl <- match(node.priors$node, names(node.states))
        for(i in 1:length(repl)) {
          if(node.priors$prior.type[i] == "normal") node.states[repl[i]] <- rnorm(1, mean = node.priors[i,2], sd = node.priors[i,3])
          if(node.priors$prior.type[i] == "uniform") node.states[repl[i]] <- runif(1, min = node.priors[i,2], max = node.priors[i,3])
          if(node.priors$prior.type[i] == "exp") node.states[repl[i]] <- node.priors[i,2] + rexp(1, rate = node.priors[i,3])
        }
      }
      d <- c(d, node.states)
    } else {
      node.states <- NULL
    }
    
    if(!is.null(root.prior)) {
      if(root.prior$prior.type == "normal") root.val <- rnorm(1, mean = root.prior[1,2], sd =root.prior[1,3])
      if(root.prior$prior.type == "uniform") root.val <- runif(1, min = root.prior[1,2], max =root.prior[1,3])
      if(root.prior$prior.type == "exp") root.val <- root.prior[1,2] + rexp(1, rate = root.prior[1,3])
    } else {
      root.val <- node.pics[1]
    }
    
    if(model == "BM" | model == "Trend") {
      inv.C <- solve(C)
    }
    
    if(model == "BM") {
      return(list("sigma" = .slater.sigsq(C, d, root.val, inv.C = inv.C)[1,1], "scaler" = NA, "root" = root.val, "node.states" = node.states, "inv.C" = inv.C))
    } else if(model == "OU") {
      return(list("sigma" = .slater.sigsq(C, d, root.val, inv.C = NULL)[1,1], "scaler" = runif(1,10^-5,1), "theta" = root.val, "root" = root.val, "node.states" = node.states))
    } else if(model == "SSP") {
      return(list("sigma" = .slater.sigsq(C, d, root.val, inv.C = NULL)[1,1], "scaler" = runif(1,10^-5,1), "root" = root.val, "node.states" = node.states))
    } else if(model == "ACDC.exp" | model == "ACDC.lin") {
      return(list("sigma" = .slater.sigsq(C, d, root.val, inv.C = NULL)[1,1], "scaler" = 0, "root" = root.val, "node.states" = node.states))
    } else if(model == "Trend") {
      return(list("sigma" = .slater.sigsq(C, d, root.val, inv.C = inv.C)[1,1], "scaler" = 0, "root" = root.val, "node.states" = node.states, "inv.C" = inv.C))
    }
  }
  
  .slater.getbrownianlk <- function(model, d, C, params) {
    if(!is.null(params$node.states)) {
      d <- c(d, params$node.states)
    }
    
    N <- length(d)
    X <- as.numeric(d[match(rownames(C), names(d))])
    root.val <- params$root
    V <- .slater.getvcv(C, model, params)
    
    if(model == "Trend") {
      EX <- root.val + (params$scaler * diag(C))
    } else if(model =="OU") {
      EX <- (root.val * exp(-params$scaler* diag(C))) + (params$theta*(1-exp(-params$scaler* diag(C))))
    } else {
      EX <- matrix(rep(root.val, N), ncol = 1)
    }
    if(model == "BM" | model == "Trend") {
      inv.V <- params$inv.C * (params$sigma ^-1)
      logL <- -t(X-EX) %*% inv.V %*% (X-EX) / 2 - N * log(2*pi) / 2 - determinant(V)$modulus[1] / 2
    } else {
      logL <- -t(X-EX) %*% solve(V) %*% (X-EX) / 2 - N * log(2*pi) / 2 - determinant(V)$modulus[1] / 2
    }
    return(as.numeric(logL))
  }
  
  .slater.initiate <- function(phy, d, C, node.priors, model, outputName, scaler.prior, root.prior, sampling.nodes) {
    params <- .slater.generatestartingvalues(phy, d, C, model, node.priors, root.prior, scaler.prior, sampling.nodes)
    
    output <- file(paste(outputName, "_model_params.txt", sep = ""), "w")
    if(model == "BM") cat("generation", "Posterior", "Prior", "logLk", "sigmasq", "root", file = output, sep = "\t")
    if(model == "OU") cat("generation", "Posterior", "Prior", "logLk", "sigmasq", "alpha", "theta", "root", file = output, sep = "\t")
    if(model == "SSP") cat("generation", "Posterior", "Prior", "logLk", "sigmasq", "alpha", "root", file = output, sep = "\t")
    if(model == "ACDC.exp") cat("generation", "Posterior", "Prior", "logLk", "sigmasq", "r", "root", file = output, sep = "\t")
    if(model == "ACDC.lin") cat("generation", "Posterior", "Prior", "logLk", "sigmasq", "beta", "root", file = output, sep = "\t")
    if(model == "Trend") cat("generation", "Posterior", "Prior", "logLk", "sigmasq", "mu", "root", file = output, sep = "\t")
    if(model == "lambda") cat("generation", "Posterior", "Prior", "logLk", "sigmasq", "lambda", "root", file = output, sep = "\t")
    cat("\n", file = output)
    
    if(sampling.nodes == TRUE) {
      node.output <- file(paste(outputName, "_nodestates.txt", sep = ""), "w")
      cat((length(phy$tip.label) + 1):(length(phy$tip.label) + phy$Nnode), file = node.output, sep = "\t")
      cat("\n", file = node.output)
    } else {
      node.output <- NULL
    }
    
    paramslogLk <- .slater.getbrownianlk(model, d, C, params)
    
    if(!is.null(node.priors)) {
      prior.prob <- .slater.dnodeprior(params$node.states[match(node.priors$node, names(params$node.states))], node.priors)
    } else {
      prior.prob <- 0
    }
    if(!is.null(root.prior)) {
      root.prior.prob <- .slater.dnodeprior(params$root, root.prior)
    } else {
      root.prior.prob <- 0
    }
    
    if(model == "OU") {
      p.out <- 4
    } else {
      p.out <- 3
    }
    
    cat(0, paramslogLk + sum(c(prior.prob, root.prior.prob)), sum(c(prior.prob, root.prior.prob)), paramslogLk, as.numeric(na.omit(as.numeric(unlist(params))[1:3])), file = output, sep = "\t")
    cat("\n", file = output)
    
    if(sampling.nodes == TRUE) {
      cat(params$root, as.numeric(na.omit(as.numeric(unlist(params$node.states)))), file = node.output, sep = "\t")
      cat("\n", file = node.output)
    }
    
    return(list(params = params, paramslogLk = paramslogLk, output = output, node.output = node.output, node.priors = prior.prob, root.prior = root.prior.prob))
  }
  
  .slater.sigsq <- function(C, d, root, inv.C = NULL) {
    N <- length(d)
    X <- as.numeric(d)
    EX <- root
    if(is.null(inv.C)) {
      inv.C <- solve(C)
    }
    return(((t(X - EX)) %*% inv.C %*% (X - EX)) / N)
  }
  
  .slater.dnodeprior <- function(fossils, node.priors) {
    pp <- numeric(nrow(node.priors))
    for(i in 1:nrow(node.priors)) {
      if(node.priors$prior.type[i] == "normal") {
        pp[i] <- dnorm(fossils[i], mean = node.priors[i,2], sd = node.priors[i,3], log = TRUE)
      }
      if(node.priors$prior.type[i] == "uniform") {
        pp[i] <- dunif(fossils[i], min = node.priors[i,2], max = node.priors[i,3], log = TRUE)
      }
      if(node.priors$prior.type[i] == "exp") {
        if(node.priors[i,2] > 0) {
          pp[i] <- dexp(fossils[i] - node.priors[i,2], rate = node.priors[i,3], log = TRUE)
        } else if(node.priors[i,2] < 0) {
          pp[i] <- dexp(fossils[i] + abs(node.priors[i,2]), rate = node.priors[i,3], log = TRUE)
        } else if(node.priors[i,2] == 0) {
          pp[i] <- dexp(fossils[i], rate = node.priors[i,3], log = TRUE)
        }
      }
    }
    return(pp)
  }
  
  .slater.updatewhichnode <- function(model, ngen, fossildata) {
    if(model == "SSP" | model == "Trend") {
      fossil.k <- 4
    } else if(model == "OU") {
      fossil.k <- 5
    } else {
      fossil.k <- 3
    }
    nfossil <- length(fossildata)
    fossil.gen <- (ngen / fossil.k)
    update.which.node <- fossil.gen %% nfossil
    if(update.which.node == 0) {
      update.which.node <- nfossil
    }
    return(update.which.node)
  }
  
  .slater.nodeupdate <- function(node, model, node.states, node.priors, pwFossil) {
    node.number <- names(node.states[node])
    if(is.null(node.number)) {
      node.number <- node.priors[1,1]
    }
    if(!node.number %in% node.priors$node) {
      return(.slater.getproposal(node.states[node], psi = pwFossil))
    } else {
      prior.row <- which(node.priors$node == node.number)
      if(node.priors$prior.type[prior.row] == "uniform") {
        return(.slater.getproposal(node.states[node], psi = pwFossil, bound = TRUE, min = node.priors[prior.row, 2], max = node.priors[prior.row, 3]))
      }
      if(node.priors$prior.type[prior.row] == "exp") {
        return(.slater.getproposal(node.states[node], psi = pwFossil, bound = TRUE, min = node.priors[prior.row, 2], max = Inf))
      }
      if(node.priors$prior.type[prior.row] == "normal") {
        return(.slater.getproposal(node.states[node], psi = pwFossil))
      }
    }
  }
  
  .slater.lkratio <- function(ProplogLk, paramslogLk) {
    return(exp(ProplogLk - paramslogLk))
  }
  
  .slater.updateacceptance <- function(k, acceptance, ngen) {
    if(k == 3) {
      if(ngen %% 3 == 1) acceptance[1] <- acceptance[1] + 1
      if(ngen %% 3 == 2) acceptance[2] <- acceptance[2] + 1
      if(ngen %% 3 == 0) acceptance[3] <- acceptance[3] + 1
    } else if (k == 4) {
      if(ngen %% 4 == 1) acceptance[1] <- acceptance[1] + 1
      if(ngen %% 4 == 2) acceptance[2] <- acceptance[2] + 1
      if(ngen %% 4 == 3) acceptance[3] <- acceptance[3] + 1
      if(ngen %% 4 == 0) acceptance[4] <- acceptance[4] + 1
    }
    return(acceptance)
  }
  
  # Preliminaries #
  td <- treedata(phy, d)
  phy <- td$phy
  d <- as.numeric(td$data)
  names(d) <- rownames(td$data)
  rm(td)
  d <- d[match(phy$tip.label, names(d))]
  
  if(is.null(node.priors)) print("no node priors have been specified. a standard model will be fit instead.")
  if (!is.null(node.priors)) {
    if(!is.data.frame(node.priors)) stop("the prior data need to be in a data frame")
    if(length(setdiff(c(as.character(node.priors[,1]), as.character(node.priors[,2])), phy$tip.label)) > 0) {
      stop("one or more node priors are defined based on taxa that are not represented in both the tree and data. Try using treedata() to check your taxon-overlap and redefine the priors")
    }
  }
  
  if(!is.null(node.priors) | sample.node.states == TRUE) {
    C <- .vcv.anc(phy, anc.nodes = TRUE)
    sampling.nodes <- TRUE
  } else {
    sampling.nodes <- FALSE
    C <- .vcv.anc(phy, anc.nodes = FALSE)
  }
  
  if(!is.null(node.priors)) {
    node.priors <- data.frame(node = apply(node.priors, 1, function(phy, x) return(.slater.mrca(phy, x[1], x[2])), phy = phy), value1 = node.priors[,3], value2 = node.priors[,4], prior.type = node.priors[,5])
    node.priors.curr.den <- numeric(nrow(node.priors))
  } else {
    node.priors.curr.den <- 0
  }
  
  if(!is.null(root.prior)) {
    if(length(root.prior) < 2) {
      stop(" if using a fossil prior on the root, you must specify distribution parameters using the root.prior argument")
    }
  }
  
  if(model == "EB") stop("EB is now an invalid option, please choose one of ACDC.lin or ACDC.exp; see documentation for more details.")
  model <- match.arg(model, c("BM", "Trend", "SSP", "ACDC.exp", "ACDC.lin", "lambda"))
  
  if(model == "ACDC.exp" | model == "ACDC.lin") {
    if(is.null(acdc.prior)) {
      scaler.prior <- .acdc.prior(phy, model)
      print(paste("no acdc.prior was specified. the values ",scaler.prior[1],":", scaler.prior[2], " have been used based on the input tree"))
    } else {
      scaler.prior <- acdc.prior
    }
  } else {
    scaler.prior <- NULL
  }
  
  if(!is.null(root.prior)) {
    root.prior <- data.frame(node = length(phy$tip.label) + 1, value1 = as.numeric(root.prior[1]), value2 = as.numeric(root.prior[2]), "prior.type" = root.prior[3])
    colnames(root.prior)[4] <- "prior.type"
  } else {
    root.prior <- NULL
  }
  
  starting <- .slater.initiate(phy, d, C, node.priors, model, outputName, scaler.prior, root.prior, sampling.nodes)
  params <- starting$params
  paramslogLk <- starting$paramslogLk
  output <- starting$output
  node.output <- starting$node.output
  curr.root.prior <- starting$root.prior
  if(!is.null(node.priors)) {
    node.priors.curr.den <- starting$node.priors
    fossil.nodes <- match(node.priors$node, names(params$node.states))
  }
  
  if(model == "BM" | model == "ACDC.lin" | model == "ACDC.exp") k <- 3
  if(model == "SSP" | model == "Trend") k <- 4
  if(model == "OU") k <- 5
  if(model == "lambda") k <- 3
  
  accept <- rep(0, k)
  
  for(ngen in 1:Ngens) {
    props <- .slater.generateproposal(ngen, propwidth, model, params, node.priors, root.prior, scaler.prior, sampling.nodes)
    ProplogLk <- .slater.getbrownianlk(model, d, C, props)
    if(!is.null(node.priors)) {
      prop.node.priors <- .slater.dnodeprior(props$node.states[fossil.nodes], node.priors)
    } else {
      prop.node.priors <- 0
    }
    if(!is.null(root.prior)) {
      prop.root.prior <- .slater.dnodeprior(props$root, root.prior)
    } else {
      prop.root.prior <- 0
    }
    
    h <- .slater.lkratio(ProplogLk, paramslogLk) * exp((prop.root.prior + sum(prop.node.priors)) - (curr.root.prior + sum(node.priors.curr.den)))
    p <- runif(1)
    
    if(h == Inf | is.na(h)) h <- 0
    if(h >= p) {
      params <- props
      paramslogLk <- ProplogLk
      curr.root.prior <- prop.root.prior
      node.priors.curr.den <- prop.node.priors
      accept <- .slater.updateacceptance(k, accept, ngen)
    }
    
    if (ngen %% sampleFreq == 0) {
      if(model == "OU") {
        p.out <- 4
      } else {
        p.out <- 3
      }
      curr.pr <- sum(c(curr.root.prior, node.priors.curr.den))
      curr.post <- paramslogLk + curr.pr
      cat(ngen, curr.post, curr.pr, paramslogLk, as.numeric(na.omit(as.numeric(unlist(params))[1:p.out])), file = output, sep = "\t")
      cat("\n", file = output)
      
      if(sampling.nodes == TRUE) {
        cat(params$root, as.numeric(na.omit(as.numeric(unlist(params$node.states)))), file = node.output, sep = "\t")
        cat("\n", file = node.output)
      }
    }
    
    if (ngen %% printFreq == 0) {
      cat("Generation:", ngen, " ", accept / (ngen / k), "\n\n")
    }
  }
  
  close(output)
  if(sampling.nodes == TRUE) {
    close(node.output)
  }
}


## simulate data under model="BM"

#x<-fastBM(tree)
#plotTree(tree)


  
  #read data
  Data <- read.csv("min20516.csv")
  View(Data)
  
  
  tree <- read.tree("min20Fixed516.nwk")
  
  
  #cut data
  
  cutData <- Data[,c(9,13),drop=FALSE] 
  cutData <- na.omit(cutData)
  cutData$Species <- gsub(" ", "_", cutData$Species) 
  includedSpecies<-cutData$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
  cutData <- cutData[!(cutData$Keep==FALSE),]
  pruned.tree$edge.length[which(pruned.tree$edge.length <=0)] <- 1e-10
  x<-cutData[,c(1,2),drop=FALSE] 
  
  all<-c(t(x))
  species<-c(t(x$Species))
  neo<-as.numeric(c(t(x$NeoplasiaPrevalence)))
  
  tree<-pruned.tree
  
  x <- setNames(neo, species)
  
  typeof(x)
  
  name.check(tree, x)


fitEB <- fitContinuous(tree, x, model="EB")
fitBM <- fitContinuous(tree, x, model="BM")
fitOU <- fitContinuous(tree, x, model="OU")
fitLambda <- fitContinuous(tree, x, model="lambda")
fitrate <- fitContinuous(tree, x, model="rate_trend")
fitmean_trend <- fitContinuous(tree, x, model="mean_trend")
fitwhite <- fitContinuous(tree, x, model="white")


evol_rate <- ace(x, tree)


#aicScores[nrow(aicScores) + 1,] = c(AIC(fitOU_sim), AIC(fitBM_sim),AIC(fitEB_sim))


#aicScores




fitOUMCMC_sim<-fitContinuousMCMC(tree, x, model="SSP", Ngens = 10000,
                  sampleFreq = 500, printFreq = 10, sample.node.states = TRUE, 
                  outputName = "ouMCMC_output")
summary(fitOUMCMC_sim)
fitBMMCMC_sim<-fitContinuousMCMC(tree, x, model="BM", Ngens = 10000,
                  sampleFreq = 500, printFreq = 10, sample.node.states = TRUE,
                  outputName = "bmMCMC_output")

fitEBMCMC_sim<-fitContinuousMCMC(tree, x, model="ACDC.exp", Ngens = 10000,
                  sampleFreq = 500, printFreq = 10,sample.node.states = FALSE,
                  outputName = "ebMCMC_output")

fitLambdaMCMC_sim<-fitContinuousMCMC(tree, x, model="lambda", Ngens = 1000,
                                 sampleFreq = 500, printFreq = 10,sample.node.states = FALSE,
                                 outputName = "lambdaMCMC_output")

OUsimtxt<-read.table("ouMCMC_output_model_params.txt", header = TRUE)


ggplot(data = OUsimtxt,aes(x= generation,y=logLk ))+
  geom_jitter()+
  scale_y_continuous(
    breaks = c(round(min(OUsimtxt$logLk)),round(median(OUsimtxt$logLk)), round(max(OUsimtxt$logLk)) ),
    labels = c(round(min(OUsimtxt$logLk)),round(median(OUsimtxt$logLk)), round(max(OUsimtxt$logLk))))+
  theme_cowplot(12)+
  geom_line()+
  labs(title = "Log Likelihood per Generation of MCMC OU")


BMsimtxt<-read.table("bmMCMC_output_model_params.txt", header = TRUE)


ggplot(data = BMsimtxt,aes(x= generation,y=logLk ))+
  geom_jitter()+
  scale_y_continuous(
    breaks = c(round(min(BMsimtxt$logLk)),round(median(BMsimtxt$logLk)), round(max(BMsimtxt$logLk)) ),
    labels = c(round(min(BMsimtxt$logLk)),round(median(BMsimtxt$logLk)), round(max(BMsimtxt$logLk))))+
  theme_cowplot(12)+
  geom_line()+
  labs(title = "Log Likelihood per Generation of MCMC BM")


EBsimtxt<-read.table("ebMCMC_output_model_params.txt", header = TRUE)


ggplot(data = EBsimtxt,aes(x= generation,y=logLk ))+
  geom_jitter()+
  theme_cowplot(12)+
  geom_line()+
  labs(title = "Log Likelihood per Generation of MCMC EB")

##this is a frankenstein of heatmap tree and labeling function should work just fine

lambdasimtxt<-read.table("lambdaMCMC_output_model_params.txt", header = TRUE)


ggplot(data = lambdasimtxt,aes(x= generation,y=logLk ))+
  geom_jitter()+
  theme_cowplot(12)+
  geom_line()+
  labs(title = "Log Likelihood per Generation of MCMC Lambda")


#subset

specData <- Data[,c(5,9,13),drop=FALSE] 
specData <- na.omit(specData)
specData$Species <- gsub(" ", "_", specData$Species) 
includedSpecies<-specData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
specData$Keep <- specData$Species %in% pruned.tree$tip.label
specData <- specData[!(specData$Keep==FALSE),]
pruned.tree$edge.length[which(pruned.tree$edge.length <=0)] <- 1e-10



mammals<-filter(specData, is.element(Clade, c("Mammalia")))
saur<-filter(specData, is.element(Clade, c("Sauropsida")))
amph<-filter(specData, is.element(Clade, c("Amphibia")))

#plot with labels



OUSimNode<-read.table("ouMCMC_output_nodestates.txt", header = TRUE)

node<-as.data.frame(OUSimNode)

nodeVal<-tail(node, n=1)

typeof(pruned.tree$node.label)
pruned.tree$edge.length


nodeVals<-c(t(nodeVal))
nodeVals<-setNames(nodeVals,seq(293,583))

typeof(nodeVals)



obj<-contMap(pruned.tree,x,plot=TRUE,method="user",anc.states=nodeVals)


## invert color scale so red is highest trait value and blue is lowest
obj<-setMap(obj,invert=TRUE)


obj
par(fg="white")

plotSimmap(obj$tree,obj$cols,type="fan", fsize = .00001 )
add.color.bar(100,cols=obj$cols,title="OU Neoplasia Prevalence", lims=range(x),digits=2)

par(fg="#3B4992FF")
arc.cladelabels(text="Mammalia",cex = .8,node=findMRCA(pruned.tree, mammals$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#EE0000FF")
arc.cladelabels(text="Sauropsida",cex = .8,node=findMRCA(pruned.tree, saur$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#008280FF")
arc.cladelabels(text="Amphibia",cex = .8,node=findMRCA(pruned.tree, amph$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)






#brownian
bmSimNode<-read.table("bmMCMC_output_nodestates.txt", header = TRUE)

node<-as.data.frame(bmSimNode)

nodeVal<-tail(node, n=1)

typeof(pruned.tree$node.label)

nodeVals<-c(t(nodeVal))
nodeVals<-setNames(nodeVals,seq(293,583))


typeof(nodeVals)



obj<-contMap(pruned.tree,x,plot=TRUE,method="user",anc.states=nodeVals)


## invert color scale so red is highest trait value and blue is lowest
obj<-setMap(obj,invert=TRUE)


obj

par(fg="white")
plotSimmap(obj$tree,obj$cols,type="fan", fsize =.00001 )
add.color.bar(100,cols=obj$cols,title="BM Neoplasia Prevalence", lims=range(x),digits=2)

par(fg="#3B4992FF")
arc.cladelabels(text="Mammalia",cex = .8,node=findMRCA(pruned.tree, mammals$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#EE0000FF")
arc.cladelabels(text="Sauropsida",cex = .8,node=findMRCA(pruned.tree, saur$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#008280FF")
arc.cladelabels(text="Amphibia",cex = .8,node=findMRCA(pruned.tree, amph$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)





