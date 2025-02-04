library(nlme)
#library(rms)
library(phytools)
library(geiger)
library(caper)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggsci)
library(patchwork)
library(poolr)
library(parallel)

library(MASS)
library(filelock)

library(R.utils)

Initialize.corPhyl=ape:::Initialize.corPhyl


vcv.corPhyl <- function(phy, data = NULL, corr = FALSE, ...) {
  if (is.null(data)) {
    labels <- attr(phy, "tree")$tip.label
    dummy.df <- data.frame(seq_along(labels), row.names = labels)
    res <- corMatrix(Initialize.corPhyl(phy, dummy.df), corr = corr)
    dimnames(res) <- list(labels, labels)
    res
  } else {
    res <- corMatrix(Initialize.corPhyl(phy, data, ...), corr = corr)
  }
  res
}

##My constants
##############
##default limits from geiger

paramDefaults <- list(
  sigsq = list(  # Variance for Brownian Motion and other models
    val = 1,
    lower = exp(-500),
    upper = exp(100)),
  
  lambda = list(  # Lambda model (Pagel 1999)
    val = 1,
    lower = exp(-500),
    upper = 1),
  
  alpha = list(  # Ornstein-Uhlenbeck (OU) model
    val = 1,
    lower = exp(-500),
    upper = exp(1)),
  
  delta = list(  # Delta model (Pagel 1999)
    val = 1,
    lower = exp(-500),
    upper = 3),
  
  kappa = list(  # Kappa model (Pagel 1999)
    val = 1,
    lower = exp(-500),
    upper = 1),
  
  a = list(  # Early burst (EB) model (rate change parameter)
    val = -0.000001,  # Default for decelerating evolution
    lower = log(10^-5),  # Based on tree depth for deceleration
    upper = 5),  # Adjust for accelerating evolution
  
  slope = list(  # Rate trend model (linear trend in rates)
    val = 0,
    lower = -100,
    upper = 100),
  
  drift = list(  # Mean trend or drift model
    val = 0,
    lower = -100,
    upper = 100)
)



##My functions
##############

var.corPhyl <- function(value, ...){
  if (!"corPhyl" %in% class(value))
    stop("The correlation object is not of class corrPhyl")
  diag(vcv.corPhyl(value, ...))
}

addVariances.phylo <- function(phy, ve){
  if (! "phylo" %in% class(phy))
    stop("phy is not of class phylo")
  if (Ntip(phy) != length(ve))
    stop("ve and phy are not compatible")
  if (!is.null(names(ve))){
    if (any(!names(ve) %in% phy$tip.label))
      stop("ve names do not match with phy names")
  } else {
    warning("ve is not named but compatible with the tree. ve will be assumed to be in phy$tip.label ordering, which is not usual")
    names(ve) <- phy$tip.label
  }
  ii <- sapply(1:Ntip(phy), function(x, e) which(e == x),
               e = phy$edge[, 2])
  phy$edge.length[ii] <- phy$edge.length[ii] + ve[phy$tip.label]
  phy
}

# Function for GLS fitting
pgls.SEy.Fixed <- function(form, data, tree, model, params = list(sigsq = 1), se = NULL, glsMethod = c("REML", "ML"), ...) {
  Call <- match.call()
  
  spp <- rownames(data)
  data <- cbind(data, spp)
  
  if (is.null(names(params)))
    stop("Params must be a named list of parameters")
  
  # Ensure sigsq and other parameter names match the expected ones
  names(params)[which(names(params) == "sig2e")] <- "sigsq"
  if (model == "rate_trend") model <- "trend"
  
  # Rescale the tree based on the selected model
  rescaledTree <- switch(model,
                         "BM" = rescale(tree, model = "BM", sigsq =params$sigsq ),
                         "lambda" = rescale(tree, model = "lambda", lambda = params$lambda, sigsq =params$sigsq),
                         "OU" = rescale(tree, model = "OU", alpha = params$alpha, sigsq =params$sigsq),
                         "delta" = rescale(tree, model = "delta", delta = params$delta, sigsq =params$sigsq),
                         "kappa" = rescale(tree, model = "kappa", kappa = params$kappa, sigsq =params$sigsq),
                         "EB" = rescale(tree, model = "EB", a = params$a, sigsq =params$sigsq),
                         "trend" = rescale(tree, model = "trend", slope = params$slope, sigsq =params$sigsq))
  rescaledTree
  if (!is.null(se)) {
    finalTree <- addVariances.phylo(rescaledTree, se^2)
  } else {
    finalTree <- rescaledTree
  }
  
  correlationStructure <- switch(model,
                                 "BM" = corBrownian(phy = finalTree, form = ~spp),
                                 "lambda" = corPagel(value = params$lambda, phy = finalTree, form = ~spp, fixed = TRUE),
                                 "OU" = corMartins(value = params$alpha, phy = finalTree, form = ~spp, fixed = TRUE),
                                 "delta" = corPagel(value = params$delta, phy = finalTree, form = ~spp, fixed = TRUE),
                                 "kappa" = corPagel(value = params$kappa, phy = finalTree, form = ~spp, fixed = TRUE),
                                 "EB" = corMartins(value = params$a, phy = finalTree, form = ~spp, fixed = TRUE),
                                 "trend" = corBrownian(phy = finalTree, form = ~spp))  # Default to corBrownian for trends
  variance <- var.corPhyl(correlationStructure, data)[spp]
  varianceStructure <- varFixed(~variance)
  data <- cbind(data, variance)
  
  # Perform GLS
  obj <- gls(form, data = data, correlation = correlationStructure, weights = varianceStructure, method = glsMethod)
  
  # Store parameters properly in the model object
  obj$phyloParams <- c(list(model = model), params)
  
  return(obj)
}

pgls.SEy.Optim <- function(form, data, model, tree, se = NULL, fixedParams = NULL, ...) {
  if (!model %in% c("BM", "lambda", "OU", "delta", "kappa", "EB", "trend"))
    stop(sprintf("pgls.SEy.Optim not implemented for %s model", model))
  if (model == "trend") model <- "rate_trend"
  
  # Extract the response variable from the formula
  response_var <- as.character(formula(form)[[2]])
  
  # Get the trait data (response variable)
  trait_data <- data[[response_var]]
  
  # Ensure the names of the trait data match the tip labels in the tree
  names(trait_data) <- rownames(data)
  
  # Set standard errors if provided, or use 0 by default
  if (is.null(se)) {
    SE <- 0
  } else {
    SE <- se
    names(SE) <- rownames(data)
  }
  
  # Fit the model without any error or warning handling
  fit <- fitContinuous(
    phy = tree,
    dat = trait_data,
    model = model,
    ...
  )
  
  # Extract estimated parameters from fitContinuous
  estimated_params <- fit$opt
  if (model == "rate_trend") model <- "trend"
  
  # Prepare parameters for pgls.SEy.Fixed based on the model
  finalParams <- switch(model,
                        "BM" = list(sigsq = fit$opt$sigsq),
                        "lambda" = list(sigsq = fit$opt$sigsq, lambda = fit$opt$lambda),
                        "OU" = list(sigsq = fit$opt$sigsq, alpha = fit$opt$alpha),
                        "delta" = list(sigsq = fit$opt$sigsq, delta = fit$opt$delta),
                        "kappa" = list(sigsq = fit$opt$sigsq, kappa = fit$opt$kappa),
                        "EB" = list(sigsq = fit$opt$sigsq, a = fit$opt$a),
                        "trend" = list(sigsq = fit$opt$sigsq, slope = fit$opt$slope))
  
  # Combine optimized parameters with any fixed parameters
  finalParams <- c(finalParams, fixedParams)
  
  # Fit the final model using the optimized parameters
  finalModel <- pgls.SEy.Fixed(form = form, data = data, tree = tree, model = model, params = finalParams, se = se)
  
  # Store the optimized parameters in the model object
  finalModel$phyloParams <- finalParams
  
  return(finalModel)
}



#make sure to run all of this before you get to work.
#pgls sey base (just run all of this)
modPgls.SEy = function (model, data, corClass = corBrownian, tree, se = NULL, 
                        method = c("REML", "ML"), interval = c(0, 1000), corClassValue=1, sig2e=NULL, ...) 
{
  Call <- match.call()
  corfunc <- corClass
  spp <- rownames(data)
  data <- cbind(data, spp)
  if (is.null(se)) 
    se <- setNames(rep(0, Ntip(tree)), tree$tip.label)[spp]
  else se <- se[spp]
  
  lk <- function(sig2e, data, tree, model, ve, corfunc, spp) {
    tree$edge.length <- tree$edge.length * sig2e
    ii <- sapply(1:Ntip(tree), function(x, e) which(e == 
                                                      x), e = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + ve[tree$tip.label]
    vf <- diag(vcv(tree))[spp]
    w <- varFixed(~vf)
    COR <- corfunc(corClassValue, tree, form = ~spp, ...)
    fit <- gls(model, data = cbind(data, vf), correlation = COR, 
               method = method, weights = w)
    -logLik(fit)
  }
  
  if (is.null(sig2e)) {
    fit <- optimize(lk, interval = interval, data = data, tree = tree, 
                    model = model, ve = se^2, corfunc = corfunc, spp = spp)
    sig2e=fit$minimum
  }
  
  tree$edge.length <- tree$edge.length * sig2e
  ii <- sapply(1:Ntip(tree), function(x, e) which(e == x), 
               e = tree$edge[, 2])
  tree$edge.length[ii] <- tree$edge.length[ii] + se[tree$tip.label]^2
  vf <- diag(vcv(tree))[spp]
  w <- varFixed(~vf)
  obj <- gls(model, data = cbind(data, vf), correlation = corfunc(corClassValue, 
                                                                  tree, form = ~spp, ...), weights = w, method = method)
  obj$call <- Call
  obj$sig2e <- sig2e
  obj
}

#Internal function
pglsSEyPagelToOptimizeLambda=function(lambda,model,data,tree,...) {
  -logLik(modPgls.SEy(model=model,data=data,tree=tree,corClassValue=lambda,corClass=corPagel,fixed=T,...)) #Returns -logLikelihood of the pgls.SEy model with lambda fixed to the value of the lambda argument. sig2e will be optimized within modPgls.SEy unless given as an argument here
}

#Function intended for users
pglsSEyPagel=function(model, data, tree, lambdaInterval=c(0,1),...){
  optimizedModel=optimize(pglsSEyPagelToOptimizeLambda,lambdaInterval,model=model,data=data,tree=tree,...) #Optimizes lambda in the lambdaInterval using the pglsSEyPagelToOptimizeLambda function
  return(modPgls.SEy(model=model,data=data,tree=tree,corClass=corPagel,fixed=T,corClassValue=optimizedModel$minimum,...)) #Returns the final model fit
}

foo <- function(x)
{
  nas <- is.na(x$coef)
  coef <- x$coef[!nas]
  cnames <- names(coef)
  coef <- matrix(rep(coef, 4), ncol = 4)
  dimnames(coef) <- list(cnames,
                         c("Estimate", "S.E.", "t", "Pr(T > |t|)"))
  df <- x$dfP - dim(coef)[1]
  coef[, 2] <- sqrt(diag(x$W))
  coef[, 3] <- coef[, 1]/coef[, 2]
  if (df < 0) {
    warning("not enough degrees of freedom to compute P-values.")
    coef[, 4] <- NA
  } else coef[, 4] <- 2 * (1 - pt(abs(coef[, 3]), df))
  coef
}

get_best_fit_model_logLik <- function(tree, neo) {
  # Ensure the names of the response data (neo) match the tip labels in the tree
  if (!all(names(neo) %in% tree$tip.label)) {
    stop("The names of the response vector (neo) must match the tip labels in the phylogenetic tree.")
  }
  
  # Define the models to test
  models_to_test <- c("BM", "lambda", "OU", "EB", "delta", "kappa", "rate_trend")
  
  # Initialize a list to store the fit results
  fit_results <- list()
  
  # Try fitting each model
  for (model in models_to_test) {
    # First attempt: fit the model normally
    fit <- fitContinuous(
      phy = tree,
      dat = neo,
      model = model
    )
    
    # Store the fit if successful
    if (!is.null(fit)) {
      fit_results[[model]] <- fit
    }
  }
  
  # Check if we have any successful fits
  if (length(fit_results) == 0) {
    stop("None of the models could be fitted successfully.")
  }
  
  # Extract log-likelihood values for each successfully fitted model
  logLik_values <- sapply(fit_results, function(fit) fit$opt$lnL, USE.NAMES = TRUE)
  
  # Get the model with the highest log-likelihood
  best_model <- names(which.max(logLik_values))
  print(best_model)
  
  return(best_model)
}

logistic_map_random_scaled <- function(iterations = 100, x0 = 0.5, r = 3.9) {
  # x0 is the initial value (between 0 and 1), r is a constant for chaos
  x <- numeric(iterations)
  x[1] <- x0
  
  # Iterate the logistic map
  for (i in 2:iterations) {
    x[i] <- r * x[i - 1] * (1 - x[i - 1])
  }
  
  # Normalize the last value to [0, 1]
  random_number <- x[iterations]
  
  # If random_number is exactly 0 or 1, adjust it slightly to avoid boundaries
  if (random_number == 0) {
    random_number <- 0.0001
  } else if (random_number == 1) {
    random_number <- 0.9999
  }
  
  # Scale it to the range [-2, 2]
  scaled_random_number <- -2 + 4 * random_number
  
  return(scaled_random_number)
}

log_error <- function(e, params, stage) {
  log_message <- sprintf("Error at stage %s with params: %s, error message: %s\n", 
                         stage, paste(params, collapse = ", "), conditionMessage(e))
  cat(log_message, file = "error_log.txt", append = TRUE)
}





dataReport <- data.frame(n = numeric(),min = numeric(), model = character(), pval = numeric(), simslope = numeric(),slope = numeric())
Data <- read.csv("min1LH.csv") # Adjust filename/path as necessary

# Define parameter settings
min_necropsies_values <- c(1,5,10,20)
#min_necropsies_values <- 20
num_species_values <- c(50)
#num_species_values <- 300
num_simulations <- 1
#slope_values <- c(0.1) 
tree <- read.tree("min1.nwk")

#slopes_check <- list()
#current_iteration = 0


no_cores <- detectCores() - 1

# Read Data
# Sample 'num_species' species from the dataset
mamData <- Data %>% 
  filter(adult_weight > 0, !is.na(adult_weight))

mamData<-filter(mamData, Class == "Mammalia")

ogData<-mamData

# Define the CSV file path
csv_file_path <- "./pglsLambdaSimResultsOU123.csv"

# Initialize the CSV file by writing the header
write.csv(data.frame(n = numeric(), min = numeric(), model = character(), pval = numeric(), simslope = numeric(), slope = numeric()), 
          file = csv_file_path, row.names = FALSE)      

run_simulation <- function(params) {
  
  
  min_necropsies <- params$min_necropsies
  num_species <- params$num_species
  seed_value <- set.seed(Sys.time() + sample(1:30000,1))  
  slope <- runif(1,min = -2, max = 2)
  sim_id <- params$sim_id
  Initialize.corPhyl=ape:::Initialize.corPhyl
  
  
  
  vcv.corPhyl <- function(phy, data = NULL, corr = FALSE, ...) {
    if (is.null(data)) {
      labels <- attr(phy, "tree")$tip.label
      dummy.df <- data.frame(seq_along(labels), row.names = labels)
      res <- corMatrix(Initialize.corPhyl(phy, dummy.df), corr = corr)
      dimnames(res) <- list(labels, labels)
      res
    } else {
      res <- corMatrix(Initialize.corPhyl(phy, data, ...), corr = corr)
    }
    res
  }
  
  ##My constants
  ##############
  ##default limits from geiger
  
  paramDefaults <- list(
    sigsq = list(  # Variance for Brownian Motion and other models
      val = 1,
      lower = exp(-500),
      upper = exp(100)),
    
    lambda = list(  # Lambda model (Pagel 1999)
      val = 1,
      lower = exp(-500),
      upper = 1),
    
    alpha = list(  # Ornstein-Uhlenbeck (OU) model
      val = 1,
      lower = exp(-500),
      upper = exp(1)),
    
    delta = list(  # Delta model (Pagel 1999)
      val = 1,
      lower = exp(-500),
      upper = 3),
    
    kappa = list(  # Kappa model (Pagel 1999)
      val = 1,
      lower = exp(-500),
      upper = 1),
    
    a = list(  # Early burst (EB) model (rate change parameter)
      val = -0.000001,  # Default for decelerating evolution
      lower = log(10^-5),  # Based on tree depth for deceleration
      upper = 5),  # Adjust for accelerating evolution
    
    slope = list(  # Rate trend model (linear trend in rates)
      val = 0,
      lower = -100,
      upper = 100),
    
    drift = list(  # Mean trend or drift model
      val = 0,
      lower = -100,
      upper = 100)
  )
  
  
  
  ##My functions
  ##############
  
  var.corPhyl <- function(value, ...){
    if (!"corPhyl" %in% class(value))
      stop("The correlation object is not of class corrPhyl")
    diag(vcv.corPhyl(value, ...))
  }
  
  addVariances.phylo <- function(phy, ve){
    if (! "phylo" %in% class(phy))
      stop("phy is not of class phylo")
    if (Ntip(phy) != length(ve))
      stop("ve and phy are not compatible")
    if (!is.null(names(ve))){
      if (any(!names(ve) %in% phy$tip.label))
        stop("ve names do not match with phy names")
    } else {
      warning("ve is not named but compatible with the tree. ve will be assumed to be in phy$tip.label ordering, which is not usual")
      names(ve) <- phy$tip.label
    }
    ii <- sapply(1:Ntip(phy), function(x, e) which(e == x),
                 e = phy$edge[, 2])
    phy$edge.length[ii] <- phy$edge.length[ii] + ve[phy$tip.label]
    phy
  }
  
  # Function for GLS fitting
  pgls.SEy.Fixed <- function(form, data, tree, model, params = list(sigsq = 1), se = NULL, glsMethod = c("REML", "ML"), ...) {
    Call <- match.call()
    
    spp <- rownames(data)
    data <- cbind(data, spp)
    
    if (is.null(names(params)))
      stop("Params must be a named list of parameters")
    
    # Ensure sigsq and other parameter names match the expected ones
    names(params)[which(names(params) == "sig2e")] <- "sigsq"
    if (model == "rate_trend") model <- "trend"
    
    # Rescale the tree based on the selected model
    rescaledTree <- switch(model,
                           "BM" = rescale(tree, model = "BM", sigsq =params$sigsq ),
                           "lambda" = rescale(tree, model = "lambda", lambda = params$lambda, sigsq =params$sigsq),
                           "OU" = rescale(tree, model = "OU", alpha = params$alpha, sigsq =params$sigsq),
                           "delta" = rescale(tree, model = "delta", delta = params$delta, sigsq =params$sigsq),
                           "kappa" = rescale(tree, model = "kappa", kappa = params$kappa, sigsq =params$sigsq),
                           "EB" = rescale(tree, model = "EB", a = params$a, sigsq =params$sigsq),
                           "trend" = rescale(tree, model = "trend", slope = params$slope, sigsq =params$sigsq))
    rescaledTree
    if (!is.null(se)) {
      finalTree <- addVariances.phylo(rescaledTree, se^2)
    } else {
      finalTree <- rescaledTree
    }
    
    correlationStructure <- switch(model,
                                   "BM" = corBrownian(phy = finalTree, form = ~spp),
                                   "lambda" = corPagel(value = params$lambda, phy = finalTree, form = ~spp, fixed = TRUE),
                                   "OU" = corMartins(value = params$alpha, phy = finalTree, form = ~spp, fixed = TRUE),
                                   "delta" = corPagel(value = params$delta, phy = finalTree, form = ~spp, fixed = TRUE),
                                   "kappa" = corPagel(value = params$kappa, phy = finalTree, form = ~spp, fixed = TRUE),
                                   "EB" = corMartins(value = params$a, phy = finalTree, form = ~spp, fixed = TRUE),
                                   "trend" = corBrownian(phy = finalTree, form = ~spp))  # Default to corBrownian for trends
    variance <- var.corPhyl(correlationStructure, data)[spp]
    varianceStructure <- varFixed(~variance)
    data <- cbind(data, variance)
    
    # Perform GLS
    obj <- gls(form, data = data, correlation = correlationStructure, weights = varianceStructure, method = glsMethod)
    
    # Store parameters properly in the model object
    obj$phyloParams <- c(list(model = model), params)
    
    return(obj)
  }
  
  pgls.SEy.Optim <- function(form, data, model, tree, se = NULL, fixedParams = NULL, ...) {
    if (!model %in% c("BM", "lambda", "OU", "delta", "kappa", "EB", "trend"))
      stop(sprintf("pgls.SEy.Optim not implemented for %s model", model))
    if (model == "trend") model <- "rate_trend"
    
    # Extract the response variable from the formula
    response_var <- as.character(formula(form)[[2]])
    
    # Get the trait data (response variable)
    trait_data <- data[[response_var]]
    
    # Ensure the names of the trait data match the tip labels in the tree
    names(trait_data) <- rownames(data)
    
    # Set standard errors if provided, or use 0 by default
    if (is.null(se)) {
      SE <- 0
    } else {
      SE <- se
      names(SE) <- rownames(data)
    }
    
    # Fit the model without any error or warning handling
    fit <- fitContinuous(
      phy = tree,
      dat = trait_data,
      model = model,
      ...
    )
    
    # Extract estimated parameters from fitContinuous
    estimated_params <- fit$opt
    if (model == "rate_trend") model <- "trend"
    
    # Prepare parameters for pgls.SEy.Fixed based on the model
    finalParams <- switch(model,
                          "BM" = list(sigsq = fit$opt$sigsq),
                          "lambda" = list(sigsq = fit$opt$sigsq, lambda = fit$opt$lambda),
                          "OU" = list(sigsq = fit$opt$sigsq, alpha = fit$opt$alpha),
                          "delta" = list(sigsq = fit$opt$sigsq, delta = fit$opt$delta),
                          "kappa" = list(sigsq = fit$opt$sigsq, kappa = fit$opt$kappa),
                          "EB" = list(sigsq = fit$opt$sigsq, a = fit$opt$a),
                          "trend" = list(sigsq = fit$opt$sigsq, slope = fit$opt$slope))
    
    # Combine optimized parameters with any fixed parameters
    finalParams <- c(finalParams, fixedParams)
    
    # Fit the final model using the optimized parameters
    finalModel <- pgls.SEy.Fixed(form = form, data = data, tree = tree, model = model, params = finalParams, se = se)
    
    # Store the optimized parameters in the model object
    finalModel$phyloParams <- finalParams
    
    return(finalModel)
  }
  
  
  
  #make sure to run all of this before you get to work.
  #pgls sey base (just run all of this)
  modPgls.SEy = function (model, data, corClass = corBrownian, tree, se = NULL, 
                          method = c("REML", "ML"), interval = c(0, 1000), corClassValue=1, sig2e=NULL, ...) 
  {
    Call <- match.call()
    corfunc <- corClass
    spp <- rownames(data)
    data <- cbind(data, spp)
    if (is.null(se)) 
      se <- setNames(rep(0, Ntip(tree)), tree$tip.label)[spp]
    else se <- se[spp]
    
    lk <- function(sig2e, data, tree, model, ve, corfunc, spp) {
      tree$edge.length <- tree$edge.length * sig2e
      ii <- sapply(1:Ntip(tree), function(x, e) which(e == 
                                                        x), e = tree$edge[, 2])
      tree$edge.length[ii] <- tree$edge.length[ii] + ve[tree$tip.label]
      vf <- diag(vcv(tree))[spp]
      w <- varFixed(~vf)
      COR <- corfunc(corClassValue, tree, form = ~spp, ...)
      fit <- gls(model, data = cbind(data, vf), correlation = COR, 
                 method = method, weights = w)
      -logLik(fit)
    }
    
    if (is.null(sig2e)) {
      fit <- optimize(lk, interval = interval, data = data, tree = tree, 
                      model = model, ve = se^2, corfunc = corfunc, spp = spp)
      sig2e=fit$minimum
    }
    
    tree$edge.length <- tree$edge.length * sig2e
    ii <- sapply(1:Ntip(tree), function(x, e) which(e == x), 
                 e = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + se[tree$tip.label]^2
    vf <- diag(vcv(tree))[spp]
    w <- varFixed(~vf)
    obj <- gls(model, data = cbind(data, vf), correlation = corfunc(corClassValue, 
                                                                    tree, form = ~spp, ...), weights = w, method = method)
    obj$call <- Call
    obj$sig2e <- sig2e
    obj
  }
  
  #Internal function
  pglsSEyPagelToOptimizeLambda=function(lambda,model,data,tree,...) {
    -logLik(modPgls.SEy(model=model,data=data,tree=tree,corClassValue=lambda,corClass=corPagel,fixed=T,...)) #Returns -logLikelihood of the pgls.SEy model with lambda fixed to the value of the lambda argument. sig2e will be optimized within modPgls.SEy unless given as an argument here
  }
  
  #Function intended for users
  pglsSEyPagel=function(model, data, tree, lambdaInterval=c(0,1),...){
    optimizedModel=optimize(pglsSEyPagelToOptimizeLambda,lambdaInterval,model=model,data=data,tree=tree,...) #Optimizes lambda in the lambdaInterval using the pglsSEyPagelToOptimizeLambda function
    return(modPgls.SEy(model=model,data=data,tree=tree,corClass=corPagel,fixed=T,corClassValue=optimizedModel$minimum,...)) #Returns the final model fit
  }
  
  foo <- function(x)
  {
    nas <- is.na(x$coef)
    coef <- x$coef[!nas]
    cnames <- names(coef)
    coef <- matrix(rep(coef, 4), ncol = 4)
    dimnames(coef) <- list(cnames,
                           c("Estimate", "S.E.", "t", "Pr(T > |t|)"))
    df <- x$dfP - dim(coef)[1]
    coef[, 2] <- sqrt(diag(x$W))
    coef[, 3] <- coef[, 1]/coef[, 2]
    if (df < 0) {
      warning("not enough degrees of freedom to compute P-values.")
      coef[, 4] <- NA
    } else coef[, 4] <- 2 * (1 - pt(abs(coef[, 3]), df))
    coef
  }
  
  get_best_fit_model_logLik <- function(tree, neo) {
    # Ensure the names of the response data (neo) match the tip labels in the tree
    if (!all(names(neo) %in% tree$tip.label)) {
      stop("The names of the response vector (neo) must match the tip labels in the phylogenetic tree.")
    }
    
    # Define the models to test
    models_to_test <- c("BM", "lambda", "OU", "EB", "delta", "kappa", "rate_trend")
    
    # Initialize a list to store the fit results
    fit_results <- list()
    
    # Try fitting each model
    for (model in models_to_test) {
      # First attempt: fit the model normally
      fit <- fitContinuous(
        phy = tree,
        dat = neo,
        model = model
      )
      
      # Store the fit if successful
      if (!is.null(fit)) {
        fit_results[[model]] <- fit
      }
    }
    
    # Check if we have any successful fits
    if (length(fit_results) == 0) {
      stop("None of the models could be fitted successfully.")
    }
    
    # Extract log-likelihood values for each successfully fitted model
    logLik_values <- sapply(fit_results, function(fit) fit$opt$lnL, USE.NAMES = TRUE)
    
    # Get the model with the highest log-likelihood
    best_model <- names(which.max(logLik_values))
    print(best_model)
    
    return(best_model)
  }
  
  # Function to write simulation results to the CSV file safely
  write_results_safely <- function(results, csv_file_path) {
    # Create a lock file to prevent race conditions
    lockfile <- paste0(csv_file_path, ".lock")
    
    # Try to acquire the lock before writing to the file
    lock <- lock(lockfile)
    on.exit(unlock(lock), add = TRUE)  # Ensure unlock happens even if there is an error
    
    # Write results in append mode
    write.table(results, file = csv_file_path, sep = ",", col.names = FALSE, append = TRUE, row.names = FALSE)
  }
  
  all_results <- list()
  compar_row <- data.frame(n = NA,min = NA, model = "Compar", pval = NA,simslope= NA,  slope = NA)
  pglsPagelOptim_row <- data.frame(n = NA,min = NA, model = "PGLSSEYOPTIM", pval = NA,simslope= NA, slope = NA)
  pglsPagel_row <- data.frame(n = NA,min = NA, model = "PGLSSEY", pval = NA,simslope= NA, slope = NA)
  
  
  
  for(sim in 1:num_simulations) {
    #current_iteration <- current_iteration + 1
    
    success <- FALSE  # Initialize a flag to check if the code runs successfully
    
    while (!success) {
      success <- tryCatch({
        # Initialize mamData
        mamData <- ogData
        
        # Sample the data while handling potential issues
        mamData <- mamData[sample(nrow(mamData), min(num_species, nrow(mamData))), ]
        n <- nrow(mamData)
        
        # Safely calculate maxNeoPrev and maxNeo
        maxNeoPrev <- max(na.omit(mamData$NeoplasiaPrevalence))
        maxNeo <- max(mamData$RecordsWithDenominators)
        
        # Safely compute mean and variance for necropsies
        mean_necropsies <- mean(mamData$Necropsies, na.rm = TRUE)
        variance_necropsies <- var(mamData$Necropsies, na.rm = TRUE)
        size_necropsy <- (mean_necropsies^2) / (variance_necropsies - mean_necropsies)
        prob_necropsy <- mean_necropsies / variance_necropsies
        
        # Define min_value and max_value
        min_value <- min_necropsies
        max_value <- maxNeo  # Ensure maxNeo is calculated correctly
        
        # Initialize storage for simulated necropsies
        simulated_necropsies <- numeric(n)
        
        # Simulate values within the specified range with error handling
        for (i in seq_len(n)) {
          max_attempts <- 1000
          attempts <- 0
          
          repeat {
            attempts <- attempts + 1
            simulated_value <- rnbinom(1, size = size_necropsy, mu = mean_necropsies)
            
            # Check if the value is within bounds
            if (simulated_value >= min_value && simulated_value <= max_value) {
              simulated_necropsies[i] <- simulated_value
              break
            }
            
            # If exceeded max attempts, use closest value and issue a warning
            if (attempts >= max_attempts) {
              cat(sprintf("Exceeded max attempts (%d) for simulation %d; using the closest valid value.\n", max_attempts, i))
              warning(sprintf("Exceeded max attempts (%d) for simulation %d; using the closest valid value.\n", max_attempts, i))
              simulated_necropsies[i] <- ifelse(simulated_value < min_value, min_value, max_value)
              break
            }
          }
        }
        
        # Assign simulated necropsies to mamData
        mamData$SimulatedNecropsies <- simulated_necropsies
        
        # Additional calculations
        mamData$SE_simpleSim <- 1 / sqrt(mamData$SimulatedNecropsies)
        
        # Handle cutData selection and preparation
        cutData <- mamData[, c(11, 9, 43, 42, 38), drop = FALSE]
        cutData$adult_weight[cutData$adult_weight < 0] <- NA
        cutData <- na.omit(cutData)
        
        # Define adult_weight
        adult_weight <- cutData$adult_weight
        adult_weightFrac <- cutData$adult_weight / max(cutData$adult_weight)
        
        # Update species names and prune the tree
        cutData$Species <- gsub(" ", "_", cutData$Species)
        includedSpecies <- cutData$Species
        pruned.tree <- drop.tip(tree, setdiff(tree$tip.label, includedSpecies))
        pruned.tree <- keep.tip(pruned.tree, pruned.tree$tip.label)
        cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
        cutData <- cutData[!(cutData$Keep == FALSE), ]
        rownames(cutData) <- cutData$Species
        # Define 'a' (ancestral state) as the lower bound but still above zero
        theta <- mean(adult_weight)
        
        
        a <- quantile(adult_weight, probs = 0.1)  # 10th percentile of adult weight
        
        # Adjust alpha and sig2 for stronger control
        alpha <- 2.5  # Higher selection strength
        sig2 <- var(adult_weight) * 0.0005  # Very low variance for minimized noise
        
        # Simulate 'x' under the OU process
        x <- fastBM(pruned.tree, a = a, theta = theta, alpha = alpha, sig2 = sig2, internal = FALSE)
        
        # Simulate 'simNeoB' under the OU model with further refined noise
        if (abs(slope - 0) < 1e-6) {
          # Minimal OU influence with lower noise
          simNeoB <- fastBM(pruned.tree, a = theta, sig2 = var(adult_weight) * 0.01, alpha = alpha)
        } else {
          # OU model with reduced noise when slope is non-zero
          noise_level <- sd(adult_weightFrac) * 0.001  # Fine-tuned noise level
          intercept <- 0.05
          
          # Add minimal noise with relationship in OU
          simNeoB <- intercept + slope * x +
            fastBM(pruned.tree, a = a, theta = theta, alpha = alpha, sig2 = sig2, internal = FALSE)
        }
        
        # Normalize and add controlled noise to simNeo values
        simNeo_normalized <- (simNeoB - min(simNeoB)) / (max(simNeoB) - min(simNeoB))
        
        # Further reduce noise level after normalization
        noise_level_normalized <- 0.001  # Minimized noise level for normalized data
        noise <- rnorm(nrow(cutData), mean = 0, sd = noise_level_normalized)*0.1
        simNeo_normalized_noisy <- simNeo_normalized + noise
        
        # Ensure values stay within the [0, 1] range
        simNeo_normalized_noisy <- pmax(0, pmin(1, simNeo_normalized_noisy))
        
        # Store values in the dataframe
        cutData$x <- x
        cutData$logx <- log10(x)
        cutData$simNeo <- simNeo_normalized
        cutData <- na.omit(cutData)
        logx <- cutData$logx
        
        # Prune tree based on available species
        pruned.tree <- drop.tip(pruned.tree, setdiff(pruned.tree$tip.label, cutData$Species))
        SE <- setNames(cutData$SE_simpleSim, cutData$Species)[rownames(cutData)]
        
        
        #Plot simNeo against x to visually inspect the relationship
        # names(simNeo_normalized)<-cutData$Species
        # fitContinuous(pruned.tree,simNeo_normalized )
        #plot(cutData$simNeo ~ cutData$logx , main = "Simulated Neoplasia (simNeo) vs x")
        # summary(pglsSEyPagel(simNeo~logx,data=cutData,tree=pruned.tree,se=SE,method = "ML"))
        # summary(pgls.SEy.Optim(simNeo~log10(x),data=cutData,tree=pruned.tree,se=SE,model = "BM"))
        # # Set success flag to TRUE if no error occurs
        TRUE
        
      }, error = function(e) {
        # In case of an error, print it and retry
        cat("Error encountered: ", conditionMessage(e), "\nRetrying...\n")
        FALSE  # Return FALSE to continue the loop and retry
      })
    }
    
    #plot(simNeo_normalized~log(x ))
    #pgls model
    tryCatch({
      p.v.adult.weight.neo<-NA
      slope.adult.weight.neo<-NA
      #cat("attempting pgls")
      adult.weight.neo<-withTimeout({pglsSEyPagel(simNeo~logx,data=cutData,tree=pruned.tree,se=SE,method = "ML")}, timeout = 500, onTimeout = "error")
      
      summary(adult.weight.neo) 
      
      #grab r squared, lambda, and p values from summary 
      
      r.v.adult.weight.neo <- summary(adult.weight.neo)$corBeta
      r.v.adult.weight.neo <- format(r.v.adult.weight.neo[2,1])
      r.v.adult.weight.neo <-signif(as.numeric(r.v.adult.weight.neo)^2, digits= 2)
      #ld.v.adult.weight.neo<- summary(adult.weight.neo)$modelStruct$corStruct
      #ld.v.adult.weight.neo <- signif(ld.v.adult.weight.neo[1], digits = 2)
      p.v.adult.weight.neo<-summary(adult.weight.neo)$tTable
      p.v.adult.weight.neo<-signif(p.v.adult.weight.neo[2,4], digits = 2)
      slope.adult.weight.neo<-summary(adult.weight.neo)$coefficients[2]
      pglsPagel_row <- data.frame(n = nrow(cutData),min = min_necropsies, model = "PGLSSEY", pval = p.v.adult.weight.neo,simslope= slope, slope = slope.adult.weight.neo)
      
    }, error = function(e) {
      # If an error occurs, assign NA to these variables
      pglsPagel_row <- data.frame(n = nrow(cutData),min = min_necropsies, model = "PGLSSEY", pval = NA, simslope = NA, slope = NA)
      
      cat("error in optim")
    })
    
    
    
    tryCatch({
      p.v.adult.weight.neoOptim<-NA
      slope.adult.weight.neoOptim<-NA
      neo<-cutData$simNeo
      names(neo)<-cutData$Species
      min_length <- 0.001  # Define the minimum branch length desired
      pruned.tree$edge.length[pruned.tree$edge.length < min_length] <- min_length
      #best_model <- get_best_fit_model_logLik(tree = pruned.tree, neo = neo)
      #print(best_model)
      
      #cat("attempting optim")
      adult.weight.neoOptim<-withTimeout({pgls.SEy.Optim(simNeo~logx,data=cutData,tree=pruned.tree,se=SE,model = "OU")}, timeout = 500, onTimeout = "error")
      
      r.v.adult.weight.neoOptim <- summary(adult.weight.neoOptim)$corBeta
      r.v.adult.weight.neoOptim <- format(r.v.adult.weight.neoOptim[2,1])
      r.v.adult.weight.neoOptim <-signif(as.numeric(r.v.adult.weight.neoOptim)^2, digits= 2)
      #ld.v.adult.weight.neoOptim<- summary(adult.weight.neoOptim)$modelStruct$corStruct
      #ld.v.adult.weight.neoOptim <- signif(ld.v.adult.weight.neoOptim[1], digits = 2)
      p.v.adult.weight.neoOptim<-summary(adult.weight.neoOptim)$tTable
      p.v.adult.weight.neoOptim<-signif(p.v.adult.weight.neoOptim[2,4], digits = 2)
      slope.adult.weight.neoOptim<-summary(adult.weight.neoOptim)$coefficients[2]
      pglsPagelOptim_row <- data.frame(n = nrow(cutData),min = min_necropsies, model = "PGLSSEYOPTIM", pval = p.v.adult.weight.neoOptim,simslope= slope, slope = slope.adult.weight.neoOptim)
    }, error = function(e) {
      # If an error occurs, assign NA to these variables
      pglsPagelOptim_row <- data.frame(n = nrow(cutData),min = min_necropsies, model = "PGLSSEYOPTIM", pval = NA,simslope= slope, slope = NA)
      cat("error in optim")
    })
    expected_neoplasia_occurrences <- cutData$simNeo * cutData$SimulatedNecropsies
    NeoplasiaOccurences <- floor(expected_neoplasia_occurrences) + 
      rbinom(nrow(cutData), 1, expected_neoplasia_occurrences - floor(expected_neoplasia_occurrences))
    # Calculate non-occurrences to ensure total matches
    NonOccurences <- cutData$SimulatedNecropsies - NeoplasiaOccurences
    #adultWeight<-log10(x)
    tryCatch({
      min_length <- 0.001  # Define the minimum branch length desired
      pruned.tree$edge.length[pruned.tree$edge.length < min_length] <- min_length
      #cat("attempting gee")
      compar <- withTimeout({suppressMessages(compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~ logx, phy = pruned.tree, family = "binomial"))}, timeout = 60, onTimeout = "error")
      
      cNpval <- foo(compar)[2,4] # Assuming summary(compar) returns the correct structure
      cNslope <- foo(compar)[2,1]
      compar_row<-data.frame(n = nrow(cutData), min = min_necropsies, model = "Compar", pval = cNpval, simslope = slope, slope = cNslope)
    }, error = function(e) {
      # If an error occurs, or if the timeout is exceeded, assign NA to these variables
      log_error(e, params, "compar.gee")
      compar_row <- data.frame(n = nrow(cutData), min = min_necropsies, model = "Compar", pval = NA, simslope = NA, slope = NA)
    })
    
    # Save the results to the CSV file after each iteration using the safe write function
    write_results_safely(rbind(pglsPagel_row, pglsPagelOptim_row, compar_row), csv_file_path)
    
    # Add the results to the in-memory list
    all_results[[length(all_results) + 1]] <- rbind(pglsPagel_row, pglsPagelOptim_row, compar_row)
  }
  return(do.call(rbind, all_results))
}



# Function to write simulation results to the CSV file safely
write_results_safely <- function(results, csv_file_path) {
  # Create a lock file to prevent race conditions
  lockfile <- paste0(csv_file_path, ".lock")
  
  # Try to acquire the lock before writing to the file
  lock <- lock(lockfile)
  on.exit(unlock(lock), add = TRUE)  # Ensure unlock happens even if there is an error
  
  # Write results in append mode
  write.table(results, file = csv_file_path, sep = ",", col.names = FALSE, append = TRUE, row.names = FALSE)
}

RNGkind("L'Ecuyer-CMRG")

set.seed(42)

cl <- makeCluster(8)

clusterSetRNGStream(cl, iseed = 42)


clusterEvalQ(cl, {
  library(nlme)
  #library(rms)
  library(phytools)
  library(geiger)
  library(caper)
  library(tidyverse)
  library(patchwork)
  library(poolr)
  library(parallel)
  library(filelock)
  library(MASS)
  library(R.utils)
})

parameter_list <- expand.grid(
  min_necropsies = min_necropsies_values,
  num_species = num_species_values
)

# Step 2: Define the number of random slopes per combination (1000 slopes)
n_slopes_per_combination <- 1000

# Step 3: Generate a vector of 1000 random slopes for each row of the parameter grid
random_slopes <- runif(nrow(parameter_list) * n_slopes_per_combination, min = -2, max = 2)

# Step 4: Replicate each combination of min_necropsies and num_species 1000 times
parameter_list <- parameter_list[rep(1:nrow(parameter_list), each = n_slopes_per_combination), ]

# Step 5: Add the random slopes to the expanded parameter list
parameter_list$random_slope <- random_slopes



clusterExport(cl, varlist = c("run_simulation", "num_simulations",
                              "min_necropsies_values", "num_species_values", 
                              "ogData", "tree",
                              "pglsSEyPagel", "pgls.SEy.Optim", "foo", "get_best_fit_model_logLik", 
                              "withTimeout","write_results_safely", "csv_file_path","parameter_list","log_error"))



# Run the simulation in parallel
results_list <- parLapply(cl, 1:nrow(parameter_list), function(i) {
  params <- parameter_list[i, ]
  run_simulation(params)
})
# Stop the cluster after the computations
stopCluster(cl)

# Combine results
final_results <- do.call(rbind, results_list)

final_results <- final_results[final_results$model %in% c("PGLSSEY", "PGLSSEYOPTIM", "Compar"), ]

# Now iterate over each column in the filtered data
for (col_name in names(final_results)) {
  # Skip the 'model' column as it is expected to be character
  if (col_name != "model") {
    # Attempt to convert the column to numeric, coercing errors to NA
    final_results[[col_name]] <- as.numeric(as.character(final_results[[col_name]]))
    
    # Check if there are any NAs introduced by coercion
    if (any(is.na(final_results[[col_name]]))) {
      message(sprintf("Non-numeric data coerced to NA in column: %s", col_name))
    }
  }
}

for (col_name in names(final_results)) {
  # Skip the 'model' column as it is expected to be character
  if (col_name != "model") {
    # Attempt to convert the column to numeric, coercing errors to NA
    final_results[[col_name]] <- as.numeric(as.character(final_results[[col_name]]))
    
    # Check if there are any NAs introduced by coercion
    if (any(is.na(final_results[[col_name]]))) {
      message(sprintf("Non-numeric data coerced to NA in column: %s", col_name))
    }
  }
}

pglsresults<-filter(final_results, model == "PGLSSEY")
pglssig<-nrow(filter(pglsresults, pval <= 0.05))
levelpgls<-pglssig/nrow(pglsresults) 

pglsOptimresults<-filter(final_results, model == "PGLSSEYOPTIM")
pglsOptimsig<-nrow(filter(pglsOptimresults, pval <= 0.05))
levelpglsOptim<-pglsOptimsig/nrow(pglsOptimresults) 


gessesults<-filter(final_results, model == "Compar")
gesssig<-nrow(filter(gessesults, pval <= 0.05))
levelgee<-gesssig/nrow(gessesults) 



write.csv(final_results, file = "./pglsLambdaSimResultsfinalOU123.csv")

