library(phytools)
library(geiger)
library(ape)
library(dplyr)

Data <- read.csv("min20516.csv")
results <- data.frame(Orders = character(), BestModel = character(), LogLik = numeric(), nSpecies = numeric())

orderCount <- table(Data$Orders)
orderCount <- orderCount[orderCount > 2]
OrderNames <- names(orderCount)

lambdasOrder<-data.frame(Order = NULL, lambda = NULL)
deltasOrder<-data.frame(Order = NULL, delta = NULL)

kappasOrder<-data.frame(Order = NULL, kappa = NULL)


for (class_name in OrderNames) {
  tree <- read.tree("min20Fixed516.nwk")
  filterData <- filter(Data, Orders == class_name)
  #for malignancy
  cutData <- filterData[, c(9,17), drop=FALSE]
  # for neo
  #cutData <- filterData[, c(9,13), drop=FALSE]
  cutData$Species <- gsub(" ", "_", cutData$Species) 
  includedSpecies <- cutData$Species
  pruned.tree <- drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree, pruned.tree$tip.label)
  cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
  cutData <- cutData[!(cutData$Keep == FALSE), ]
  rownames(cutData) <- cutData$Species
  pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10
  x <- setNames(as.numeric(cutData[, 2]), cutData[, 1])
  
  name.check(pruned.tree, x)
  
  results <- tryCatch({
    fitEB_sim <- fitContinuous(pruned.tree, x, model = "EB")
    fitBM_sim <- fitContinuous(pruned.tree, x, model = "BM")
    fitOU_sim <- fitContinuous(pruned.tree, x, model = "OU")
    fitRateTrend_sim <- fitContinuous(pruned.tree, x, model = "rate_trend")
    fitMeanTrend_sim <- fitContinuous(pruned.tree, x, model = "mean_trend")
    fitWhite_sim <- fitContinuous(pruned.tree, x, model = "white")
    fitLambda_sim <- fitContinuous(pruned.tree, x, model="lambda")
    
    fitKappa_sim <- fitContinuous(pruned.tree, x, model="kappa")
    
    fitDelta_sim <- fitContinuous(pruned.tree, x, model="delta")
    
    lambdarow<-data.frame(group = class_name, lambda = fitLambda_sim$opt$lambda)
    lambdasOrder<-rbind(lambdasOrder,lambdarow)
    
    deltarow<-data.frame(group = class_name, delta = fitDelta_sim$opt$delta)
    deltasOrder<-rbind(deltasOrder,deltarow)
    
    kapparow<-data.frame(group = class_name, kappa = fitKappa_sim$opt$kappa)
    kappasOrder<-rbind(kappasOrder,kapparow)
    
    
    testResults <- rbind(
      data.frame(Orders = class_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "white", LogLik = fitWhite_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "lambda", LogLik = fitLambda_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "kappa", LogLik = fitKappa_sim$opt$lnL, nSpecies = length(x)),
      data.frame(Orders = class_name, BestModel = "delta", LogLik = fitDelta_sim$opt$lnL, nSpecies = length(x))
    )
    max_row <- testResults[which.max(testResults$LogLik), ]
    rbind(results, max_row)
  }, error = function(e) {
    cat("Error encountered: ", e$message, "\nAttempting to reorder tree and retry.\n")
    reordered_tree <- reorder.phylo(pruned.tree, order = "postorder")
    tryCatch({
      fitEB_sim <- fitContinuous(reordered_tree, x, model = "EB")
      fitBM_sim <- fitContinuous(reordered_tree, x, model = "BM")
      fitOU_sim <- fitContinuous(reordered_tree, x, model = "OU")
      fitRateTrend_sim <- fitContinuous(reordered_tree, x, model = "rate_trend")
      fitMeanTrend_sim <- fitContinuous(reordered_tree, x, model = "mean_trend")
      fitWhite_sim <- fitContinuous(reordered_tree, x, model = "white")
      fitLambda_sim <- fitContinuous(reordered_tree, x, model="lambda")
      fitKappa_sim <- fitContinuous(reordered_tree, x, model="kappa")
      fitDelta_sim <- fitContinuous(reordered_tree, x, model="delta")
      
      
      
      testResults <- rbind(
        data.frame(Orders = class_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "white", LogLik = fitWhite_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "lambda", LogLik = fitLambda_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "kappa", LogLik = fitKappa_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Orders = class_name, BestModel = "delta", LogLik = fitDelta_sim$opt$lnL, nSpecies = length(x))
        
      )
      max_row <- testResults[which.max(testResults$LogLik), ]
      rbind(results, max_row)
    }, error = function(e) {
      cat("Second attempt failed: ", e$message, "\nReturning NA for all result rows.\n")
      NAresult <- data.frame(Orders = class_name, BestModel = NA, LogLik = NA, nSpecies = length(x))
      rbind(results, NAresult)
    })
  })
}

# Print the final results to check
print(results)

modelsBestFit<-na.omit(results)






# Filtering the data
filtered_data <- Data %>% 
  filter(Orders %in% as.list(modelsBestFit$Orders))

bestFitSpecies<-left_join(filtered_data, modelsBestFit, by = "Orders")

cutDataOrder<-bestFitSpecies[c(9,6,44)]



cutDataOrder<-bestFitSpecies[c(9,6,44)]
tree <- read.tree("min20Fixed516.nwk")

cutDataOrder$Species <- gsub(" ", "_", cutDataOrder$Species) 
includedSpecies<-cutDataOrder$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutDataOrder$Keep <- cutDataOrder$Species %in% pruned.tree$tip.label
cutDataOrder <- cutDataOrder[!(cutDataOrder$Keep==FALSE),]
rownames(cutDataOrder)<-cutDataOrder$Species


speciesModelsOrder <- setNames((cutDataOrder[, 3]), (cutDataOrder[, 1]))

ou <- names(speciesModelsOrder)[speciesModelsOrder == "OU"]
EB <- names(speciesModelsOrder)[speciesModelsOrder == "EB"]
BM <- names(speciesModelsOrder)[speciesModelsOrder == "BM"]
white <- names(speciesModelsOrder)[speciesModelsOrder == "white"]
rate_trend <- names(speciesModelsOrder)[speciesModelsOrder == "rate_trend"]
lambda <- names(speciesModelsOrder)[speciesModelsOrder == "lambda"]






#for Class now 



Data <- read.csv("min20516.csv")
resultsClass <- data.frame(Class = character(), BestModel = character(), LogLik = numeric(), nSpecies = numeric())

ClassCount <- table(Data$Class)
ClassCount <- ClassCount[ClassCount > 2]
ClassNames <- names(ClassCount)

lambdasClass<-data.frame(Class = NULL, lambda = NULL)
kappasClass<-data.frame(Class = NULL, kappa = NULL)
deltasClass<-data.frame(Class = NULL, delta = NULL)


for (class_name in ClassNames) {
  tree <- read.tree("min20Fixed516.nwk")
  filterData <- filter(Data, Class == class_name)
  #for malignancy
  cutData <- filterData[, c(9,17), drop=FALSE]
  # for neo
  #cutData <- filterData[, c(9,13), drop=FALSE]
  cutData$Species <- gsub(" ", "_", cutData$Species) 
  includedSpecies <- cutData$Species
  pruned.tree <- drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree, pruned.tree$tip.label)
  cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
  cutData <- cutData[!(cutData$Keep == FALSE), ]
  rownames(cutData) <- cutData$Species
  pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10
  x <- setNames(as.numeric(cutData[, 2]), cutData[, 1])
  
  name.check(pruned.tree, x)
  
  tryCatch({
      fitEB_sim <- fitContinuous(pruned.tree, x, model = "EB")
      fitBM_sim <- fitContinuous(pruned.tree, x, model = "BM")
      fitOU_sim <- fitContinuous(pruned.tree, x, model = "OU")
      fitLambda_sim <- fitContinuous(pruned.tree, x, model="lambda")
      fitKappa_sim <- fitContinuous(pruned.tree, x, model="kappa")
      fitDelta_sim <- fitContinuous(pruned.tree, x, model="delta")
      fitRateTrend_sim <- fitContinuous(pruned.tree, x, model = "rate_trend")
      fitMeanTrend_sim <- fitContinuous(pruned.tree, x, model = "mean_trend")
      fitWhite_sim <- fitContinuous(pruned.tree, x, model = "white")
      lambdarow<-data.frame(Class = class_name, lambda = fitLambda_sim$opt$lambda)
      lambdasClass<-rbind(lambdasClass,lambdarow)
      
      kapparow<-data.frame(Class = class_name, kappa = fitKappa_sim$opt$kappa)
      kappasClass<-rbind(kappasClass,kapparow)
      
      deltarow<-data.frame(Class = class_name, delta = fitDelta_sim$opt$delta)
      deltasClass<-rbind(deltasClass,deltarow)
      
      
      testResults <- rbind(
        data.frame(Class = class_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "white", LogLik = fitWhite_sim$opt$lnL, nSpecies = length(x)),        
        data.frame(Class = class_name, BestModel = "lambda", LogLik = fitLambda_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "kappa", LogLik = fitKappa_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "delta", LogLik = fitDelta_sim$opt$lnL, nSpecies = length(x))
      )
    max_row <- testResults[which.max(testResults$LogLik), ]
    resultsClass<-rbind(resultsClass, max_row)
  }, error = function(e) {
    cat("Error encountered: ", e$message, "\nAttempting to reorder tree and retry.\n")
    reordered_tree <- reorder.phylo(pruned.tree, order = "postorder")
    tryCatch({
      fitEB_sim <- fitContinuous(reordered_tree, x, model = "EB")
      fitBM_sim <- fitContinuous(reordered_tree, x, model = "BM")
      fitOU_sim <- fitContinuous(reordered_tree, x, model = "OU")
      fitLambda_sim <- fitContinuous(pruned.tree, x, model="lambda")
      fitKappa_sim <- fitContinuous(pruned.tree, x, model="kappa")
      fitDelta_sim <- fitContinuous(pruned.tree, x, model="delta")
      fitRateTrend_sim <- fitContinuous(reordered_tree, x, model = "rate_trend")
      fitMeanTrend_sim <- fitContinuous(reordered_tree, x, model = "mean_trend")
      fitWhite_sim <- fitContinuous(reordered_tree, x, model = "white")

      
      
      testResults <- rbind(
        data.frame(Class = class_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "white", LogLik = fitWhite_sim$opt$lnL, nSpecies = length(x)),        
        data.frame(Class = class_name, BestModel = "EB", LogLik = fitLambda_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "lambda", LogLik = fitLambda_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "kappa", LogLik = fitKappa_sim$opt$lnL, nSpecies = length(x)),
        data.frame(Class = class_name, BestModel = "delta", LogLik = fitDelta_sim$opt$lnL, nSpecies = length(x))
        
      )
      max_row <- testResults[which.max(testResults$LogLik), ]
      resultsClass<-rbind(resultsClass, max_row)
    }, error = function(e) {
      cat("Second attempt failed: ", e$message, "\nReturning NA for all result rows.\n")
      NAresult <- data.frame(Class = class_name, BestModel = NA, LogLik = NA, nSpecies = length(x))
      resultsClass<-rbind(resultsClass, NAresult)
    })
  })
}


modelsBestFitClass<-na.omit(resultsClass)

# Filtering the data
filtered_data <- Data %>% 
  filter(Class %in% as.list(modelsBestFitClass$Class))

bestFitSpeciesClass<-left_join(filtered_data, modelsBestFitClass, by = "Class")

cutDataClass<-bestFitSpeciesClass[c(9,4,44)]
cutDataClass$Species <- gsub(" ", "_", cutDataClass$Species) 

# Get unique species from DataFrame A
# Get unique species from cutDataOrder
unique_species_order <- unique(cutDataOrder$Species)

# Filter cutDataClass to only include rows with species present in cutDataOrder
cutDataClass<- cutDataClass[cutDataClass$Species %in% unique_species_order, ]


tree <- read.tree("min20Fixed516.nwk")

#cutDataClass$Species <- gsub(" ", "_", cutDataClass$Species) 
includedSpecies<-cutDataClass$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutDataClass$Keep <- cutDataClass$Species %in% pruned.tree$tip.label
cutDataClass <- cutDataClass[!(cutDataClass$Keep==FALSE),]
rownames(cutDataClass)<-cutDataClass$Species


speciesModelsClass <- setNames((cutDataClass[, 3]), (cutDataClass[, 1]))


ou <- names(speciesModelsClass)[speciesModelsClass == "OU"]
EB <- names(speciesModelsClass)[speciesModelsClass == "EB"]
BM <- names(speciesModelsClass)[speciesModelsClass == "BM"]
white <- names(speciesModelsClass)[speciesModelsClass == "white"]
rate_trend <- names(speciesModelsClass)[speciesModelsClass == "rate_trend"]
lambda <- names(speciesModelsClass)[speciesModelsClass == "lambda"]


Data <- read.csv("min20516.csv")

tree <- read.tree("min20Fixed516.nwk")
#for malignancy
cutDataAll <- Data[, c(9,17, 6), drop=FALSE]
cutDataAll$Species <- gsub(" ", "_", cutDataAll$Species) 

cutDataAll<- cutDataAll[cutDataAll$Species %in% unique_species_order, ]



# for neo
#cutData <- filterData[, c(9,13), drop=FALSE]

includedSpecies <- cutDataAll$Species
pruned.tree <- drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree, pruned.tree$tip.label)
cutDataAll$Keep <- cutDataAll$Species %in% pruned.tree$tip.label
cutDataAll <- cutDataAll[!(cutDataAll$Keep == FALSE), ]
rownames(cutDataAll) <- cutDataAll$Species
pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10
x <- setNames(as.numeric(cutDataAll[, 2]), cutDataAll[, 1])

name.check(pruned.tree, x)

lambdasAll<-data.frame(Class = NULL, lambda = NULL)
kappasAll<-data.frame(Class = NULL, kappa = NULL)
deltasAll<-data.frame(Class = NULL, delta = NULL)

  fitEB_simAll <- fitContinuous(pruned.tree, x, model = "EB")
  fitBM_simAll <- fitContinuous(pruned.tree, x, model = "BM")
  fitOU_simAll <- fitContinuous(pruned.tree, x, model = "OU")
  fitLambdaAll<-fitContinuous(pruned.tree, x, model = "lambda")
  lambdarow<-data.frame(group = "all", lambda = fitLambdaAll$opt$lambda)
  lambdasAll<-rbind(lambdasAll,lambdarow)
  fitKappaAll<-fitContinuous(pruned.tree, x, model = "kappa")
  kapparow<-data.frame(group = "all", kappa = fitKappaAll$opt$kappa)
  kappasAll<-rbind(kappasAll,kapparow)
  fitDeltaAll<-fitContinuous(pruned.tree, x, model = "delta")
  deltarow<-data.frame(group = "all", delta = fitDeltaAll$opt$delta)
  deltasAll<-rbind(deltasAll,deltarow)
  fitRateTrend_simAll <- fitContinuous(pruned.tree, x, model = "rate_trend")
  fitMeanTrend_simAll <- fitContinuous(pruned.tree, x, model = "mean_trend")
  fitWhite_simAll <- fitContinuous(pruned.tree, x, model = "white")
  
  testResultsAll <- rbind(
    data.frame(BestModel = "EB", LogLik = fitEB_simAll$opt$lnL, nSpecies = length(x)),
    data.frame( BestModel = "BM", LogLik = fitBM_simAll$opt$lnL, nSpecies = length(x)),
    data.frame( BestModel = "OU", LogLik = fitOU_simAll$opt$lnL, nSpecies = length(x)),
    data.frame( BestModel = "rate_trend", LogLik = fitRateTrend_simAll$opt$lnL, nSpecies = length(x)),
    data.frame(BestModel = "mean_trend", LogLik = fitMeanTrend_simAll$opt$lnL, nSpecies = length(x)),
    data.frame(BestModel = "white", LogLik = fitWhite_simAll$opt$lnL, nSpecies = length(x)),
    data.frame(BestModel = "lambda", LogLik = fitLambdaAll$opt$lnL, nSpecies = length(x)),
    data.frame(BestModel = "kappa", LogLik = fitKappaAll$opt$lnL, nSpecies = length(x)),
    data.frame(BestModel = "delta", LogLik = fitDeltaAll$opt$lnL, nSpecies = length(x))
    
  )
    max_row <- testResultsAll[which.max(testResultsAll$LogLik), ]
  resultsAll<-rbind(testResultsAll, max_row)


  modelsBestFitAll<-na.omit(resultsAll)
  
  


name.check(pruned.tree, speciesModelsClass)
reordered_tree <- reorder.phylo(pruned.tree, order = "cladewise")
painted <- reordered_tree

#painted <- paintSubTree(painted, node=findMRCA(painted, cutDataAll$Species), state="OU", anc="0")

painted <- paintSubTree(painted, node=findMRCA(painted, cutDataAll$Species), state="Pagel's Lambda", anc="0")


# Loop through each row in the dataframe
for (i in 1:nrow(modelsBestFitClass)) {
  Class_data <- filter(cutDataClass, Class == modelsBestFitClass$Class[i])
  
  # Check the best model type and apply the corresponding functions
  if (modelsBestFitClass$BestModel[i] == "OU") {
    if (Class_data$Class[i] == "Reptilia"){
      painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species)+1, state="OU", anc="0")}
    else{
      painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="OU", anc="0")}
  } else if (modelsBestFitClass$BestModel[i] == "EB") {
    painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="EB", anc="0")
  } else if (modelsBestFitClass$BestModel[i] == "BM") {
    painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="BM", anc="0")
  }else if (modelsBestFitClass$BestModel[i] == "white") {
    painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="white", anc="0")
  }else if (modelsBestFitClass$BestModel[i] == "rate_trend") {
    painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="Rate Trend", anc="0")
  }else if (modelsBestFitClass$BestModel[i] == "lambda") {
    classLambda<-left_join(Class_data, lambdasClass, by = "Class")
      #painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state=paste("Lambda~",round(classLambda$lambda[1],1) ), anc="0")
      painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="Pagel's Lambda", anc="0")
  }else if (modelsBestFitClass$BestModel[i] == "kappa") {
    classKappa<-left_join(Class_data, kappasClass, by = "Class")
    painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="Pagel's Kappa", anc="0")
  }else if (modelsBestFitClass$BestModel[i] == "delta") {
    classDelta<-left_join(Class_data, deltasClass, by = "Class")
    painted <- paintSubTree(painted, node=findMRCA(painted, Class_data$Species), state="Pagel's Delta", anc="0")
  }
  }
    
  



# Loop through each row in the dataframe
for (i in 1:nrow(modelsBestFit)) {
  order_data <- filter(cutDataOrder, Orders == modelsBestFit$Orders[i])
  colnames(order_data)<-c("Species", "group", "BestModel","Keep")
  
  # Check the best model type and apply the corresponding functions
  if (modelsBestFit$BestModel[i] == "OU") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="OU", anc="0")
  } else if (modelsBestFit$BestModel[i] == "EB") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="EB", anc="0")
  } else if (modelsBestFit$BestModel[i] == "BM") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="BM", anc="0")
  }else if (modelsBestFit$BestModel[i] == "white") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="white", anc="0")
  }else if (modelsBestFit$BestModel[i] == "rate_trend") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="Rate Trend", anc="0")
  }else if (modelsBestFit$BestModel[i] == "lambda") {
    orderLambda<-left_join(order_data, lambdasOrder, by = "group")
    #painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state=paste("Pagel's Lambda~",round(orderLambda$lambda[1],1) ), anc="0")
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="Pagel's Lambda", anc="0")
  }else if (modelsBestFit$BestModel[i] == "delta") {
    orderDelta<-left_join(order_data, deltasOrder, by = "group")
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="Pagel's Delta", anc="0")
  }else if (modelsBestFit$BestModel[i] == "kappa") {
    orderKappa<-left_join(order_data, kappasOrder, by = "group")
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="Pagel's Kappa", anc="0")
  }
}

names(painted$maps[1])
first_element_names <- lapply(painted$maps, function(x) names(x)[1])

# Printing the names of the first elements
unique(first_element_names)


colors <- c("black", "#CC5500", "#008080", "gold", "#2E4C6D", "pink")
names(colors)<-c("0",unique(first_element_names))
plotSimmap(painted, lwd=2, split.vertical=TRUE, ftype="i", fsize = 0.0000001, colors = colors)




# Loop through each row in the dataframe
for (i in 1:nrow(modelsBestFit)) {
  order_data <- filter(cutDataAll, Orders == modelsBestFit$Orders[i])
  cladelabels(painted, text = modelsBestFit$Orders[i], node = findMRCA(painted, order_data$Species), offset = 1, orientation = "horizontal")
}


add.simmap.legend(leg=names(colors)[2:6], colors=colors[2:6])
