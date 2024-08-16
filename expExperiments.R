library(MASS)
library(nlme)
#library(rms)
library(phytools)
library(geiger)
library(caper)
library(tidyverse)
library(parallel)
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
a0 <- exp(-500)
#aPosI <- exp(100) #Such a high positive infinite limit was generating issues, probably during tree transformation
aPosI <- exp(10)

#Parameter defaults
paramDefaults <- list(
  sigsq = list(
    val = 1,
    lower = a0,
    upper = aPosI),
  lambda = list(
    val = 1,
    lower = a0,
    upper = 1),
  alpha = list(
    val = 1,
    lower = a0,
    upper = exp(1))
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

pgls.SEy.Fixed <- function(form, data, tree, model, params = list(sigsq = 1), se = NULL, glsMethod = c("REML", "REML"), ...) { ##TODO double-check the ellipsis. Do I need it?
  Call <- match.call()
  
  ##This function works with ordering following rownames(data)
  spp <- rownames(data)
  data <- cbind(data, spp)
  
  if(is.null(names(params)))
    stop("Params must be a named list of parameters")
  
  #sig2e and sigsq are synonyms, only using one
  names(params)[which(names(params)=="sig2e")]="sigsq"
  rescaledTree <- do.call(rescale, c(list(x=tree,model=model),params))
  
  #Add VE to tree
  if(!is.null(se)){
    finalTree <- addVariances.phylo(rescaledTree,se^2)
  } else {
    finalTree <- rescaledTree
  }
  
  #Use tree to generate FIXED COR and VAR
  correlationStructure=corBrownian(phy = finalTree, form = ~spp)
  variance=var.corPhyl(correlationStructure,data)[spp]
  varianceStructure=varFixed(~variance)
  data=cbind(data,variance)
  
  #gls
  obj <- gls(form, data = data, correlation = correlationStructure, weights = varianceStructure, method = glsMethod)
  #TODO This is generating ugly outputs. May want to work on it in the future
  #obj$call <- Call 
  obj$phyloParams <- c(list(model = model),params)
  obj
}


pgls.SEy.Optim <- function(form, data, model, tree, se = NULL, fixedParams = NULL, ...) {
  if(!model %in% c("BM","lambda","OU"))
    stop(sprintf("pgls.SEy.Optim not implemented for %s model",model))
  
  dataArgs <- list(
    form = form,
    data = data,
    model = model,
    tree = tree,
    se = se
  )
  
  lk <- function(parVals, optimParams, fixedParams = NULL, ...) {
    if(length(optimParams) != length(parVals))
      stop("Number of parameters not compatible with the model specified")
    params <- as.list(parVals)
    names(params) <- optimParams
    -logLik(pgls.SEy.Fixed(..., params = c(params,fixedParams)))
  } 
  
  if(model == "BM"){
    optimParams <- c("sigsq")
  } else if (model == "lambda" | model == "Pagel") {
    dataArgs$model <- "lambda"
    optimParams <- c("sigsq", "lambda")
  } else if (model == "OU" | model == "alpha") {
    dataArgs$model <- "OU"
    optimParams <- c("sigsq", "alpha")
  }
  
  for (fixedParam in names(fixedParams)){
    optimParams <- optimParams[-which(optimParams==fixedParam)]
  }
  
  if(is.null(se) & is.null(fixedParams[["sigsq"]])) {
    warning("Sigsq can't be estimated if SE is not provided")
    optimParams <- optimParams[-which(optimParams=="sigsq")]
    fixedParams <- list(sigsq = 1)
  }
  
  #Setting up optimization parameters depending on the model selected and provided data
  optimArgs <- list(
    par = NULL,
    optimParams = optimParams,
    fixedParams = fixedParams,
    fn = lk,
    lower = NULL,
    upper = NULL
  )
  
  for (parameter in optimParams){
    thisParamDefaults <- paramDefaults[[parameter]]
    if(is.null(thisParamDefaults))
      stop(sprintf("Parameter %s not recognized",parameter))
    
    optimArgs$par <- c(optimArgs$par,thisParamDefaults$val)
    optimArgs$lower <- c(optimArgs$lower,thisParamDefaults$lower)
    optimArgs$upper <- c(optimArgs$upper,thisParamDefaults$upper)
  }
  
  bmArgs <- list()
  
  optimlk <- NULL
  if(length(optimArgs$optimParams)>0) {
    if(length(optimArgs$optimParams)>1) {
      optimArgs$method <- "L-BFGS-B"
    } else {
      optimArgs$method <- "Brent"
    }
    
    optimResults <- do.call(optim,
                            c(optimArgs, dataArgs, ...))
    
    params <- as.list(optimResults$par)
    names(params) <- optimParams
    
    bmArgs <- c(list(params = c(params,fixedParams)))
    optimlk <- optimResults$value
    
  }
  finalModel <- do.call(pgls.SEy.Fixed,
                        c(dataArgs, bmArgs, ...))
  finallk <- as.numeric(logLik(finalModel))
  
  if ((!is.null(optimlk)) && abs(finallk+optimlk)>1E-6)
    stop(sprintf("likelihood issue, optimization %f, final %f",optimlk, finallk))
  
  finalModel
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


dataReport <- data.frame(n = numeric(),min = numeric(), model = character(), pval = numeric(), simslope = numeric(),slope = numeric())
Data <- read.csv("min1LH.csv") # Adjust filename/path as necessary

# Define parameter settings
min_necropsies_values <- c(1,5,10,20)
#min_necropsies_values <- 20
num_species_values <- c(50, 200)
#num_species_values <- 300
num_simulations <- 1
tree <- read.tree("min1.nwk")


slope_values <- c(.021, 0, -.2)  # Array of slope values


no_cores <- detectCores() - 1

# Read Data
      # Sample 'num_species' species from the dataset
      mamData <- Data %>% 
        filter(adult_weight > 0, !is.na(adult_weight))
      
      
      total_iterations <- num_simulations * length(slope_values) * length(num_species_values) * length(min_necropsies_values)
      current_iteration <- 0
      
      # Initialize progress tracking
      progress_threshold <- 10
      next_progress_update <- progress_threshold
      

run_simulation <- function(params) {
  compar_row <- data.frame(n = NA,min = NA, model = "Compar", pval = NA,simslope= NA,  slope = NA)
  # Extract the slope directly from params as it's directly passed
  slope <- params
  
  # Initialize an empty list to store results of each simulation
  all_results <- list()
  
  
  for(min_necropsies in min_necropsies_values) {
    for(num_species in num_species_values) {
      for(sim in 1:num_simulations) {
        current_iteration <- current_iteration + 1
        progress <- (current_iteration / total_iterations) * 100
        
        # Check if progress has reached the next 10% threshold
        if (progress >= next_progress_update) {
          cat(sprintf("Progress: [%-10s] %d%%\n", paste(rep("=", next_progress_update / 10), collapse = ""), next_progress_update))
          next_progress_update <<- next_progress_update + progress_threshold
        }

        
      mamData <- mamData[sample(nrow(mamData), min(num_species, nrow(mamData))), ]
      n <- nrow(mamData)
      maxNeoPrev<-max(na.omit(mamData$NeoplasiaPrevalence))
      maxNeo<-max(mamData$NeoplasiaWithDenominators)
      
      mean_necropsies <- mean(mamData$Necropsies, na.rm = TRUE)
      variance_necropsies <- var(mamData$Necropsies, na.rm = TRUE)
      size_necropsy <- (mean_necropsies^2) / (variance_necropsies - mean_necropsies)
      prob_necropsy <- mean_necropsies / variance_necropsies
      
      # Specify your minimum and maximum constraints
      min_value <- min_necropsies
      max_value <- maxNeo  # Ensure you've defined or calculated maxNeo appropriately
      
      # Initialize the storage for simulated values
      simulated_necropsies <- numeric(n)  # n should be defined as the number of simulations you want
      
      # Simulate values within the specified range
      for (i in seq_len(n)) {
        repeat {
          simulated_value <- rnbinom(1, size = size_necropsy, mu = mean_necropsies)
          if (simulated_value >= min_value && simulated_value <= max_value) {
            simulated_necropsies[i] <- simulated_value
            break
          }
        }
      }
      
      mamData$SimulatedNecropsies <- simulated_necropsies



#adult weight models
#adult weight neo

mamData$SE_simpleSim <-1/sqrt(mamData$SimulatedNecropsies)

cutData <- mamData[,c(11,9,43,42,38),drop=FALSE] 
cutData$adult_weight[cutData$adult_weight < 0] <- NA
cutData <- na.omit(cutData)


cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simpleSim,cutData$Species)[rownames(cutData)]

# Adjusting the simulation of 'simNeo' with increased variability
adult_weightFrac<-cutData$adult_weight/max(cutData$adult_weight)
noise_level <- sd(adult_weightFrac)*0.0000001 # Keep noise level consistent
intercept<-min(adult_weightFrac)
simNeo <- intercept+slope * adult_weightFrac + rnorm(nrow(cutData), mean = mean(adult_weightFrac), sd = noise_level)
simNeo[simNeo < 0] <- 0
cutData$simNeo<-simNeo
#pgls model
adult.weight.neo<-pglsSEyPagel(simNeo~log10(adult_weight),data=cutData,tree=pruned.tree,se=SE,method = "ML")

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




tryCatch({
  adult.weight.neoOptim<-pgls.SEy.Optim(simNeo~log10(adult_weight),data=cutData,tree=pruned.tree,se=SE,glsMethod = "ML",model = "OU")
  
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


NeoplasiaOccurences<-round(as.numeric(cutData$simNeo*cutData$SimulatedNecropsies))
NonOccurences<-round(cutData$SimulatedNecropsies- NeoplasiaOccurences)
adultWeight<-log10(cutData$adult_weight)
pglsPagel_row <- data.frame(n = nrow(cutData),min = min_necropsies, model = "PGLSSEY", pval = p.v.adult.weight.neo,simslope= slope, slope = slope.adult.weight.neo)
tryCatch({
  compar_row <- withTimeout({
    compar <- suppressMessages(compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~ adultWeight, phy = pruned.tree, family = "binomial"))
    cNpval <- foo(compar)[2,4] # Assuming summary(compar) returns the correct structure
    cNslope <- foo(compar)[2,1]
    data.frame(n = nrow(cutData), min = min_necropsies, model = "Compar", pval = cNpval, simslope = slope, slope = cNslope)
  }, timeout = 60, onTimeout = "error")
}, error = function(e) {
  # If an error occurs, or if the timeout is exceeded, assign NA to these variables
  cat("gee time out/ error")
  compar_row <- data.frame(n = nrow(cutData), min = min_necropsies, model = "Compar", pval = NA, simslope = NA, slope = NA)
})
# ggplot(cutData, aes(x = adultWeight, y = NeoplasiaOccurences / (NeoplasiaOccurences + NonOccurences))) +
#   geom_point() +
#   labs(title = "Effect of Adult Weight on Neoplasia Occurrences", x = "Adult Weight", y = "Probability of Neoplasia Occurrence") +
#   theme_minimal()


all_results[[length(all_results) + 1]] <- rbind(pglsPagel_row,pglsPagelOptim_row, compar_row)

# Combine results for this simulation
#dataReport <- rbind(dataReport, pglsPagel_row, compar_row)

      }
    }
  }
  result_rows <- do.call(rbind, all_results)
  return(result_rows)
}




results_list <- mclapply(slope_values, run_simulation, mc.cores = no_cores)

# Combine all results from each slope into one data frame
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



write.csv(final_results, file = "./pglsOUSimResults812.csv")

