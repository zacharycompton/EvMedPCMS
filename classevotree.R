library(phytools)
library(geiger)
library(ape)
library(dplyr)

Data <- read.csv("min20516.csv")
results <- data.frame(Orders = character(), BestModel = character(), LogLik = numeric(), nSpecies = numeric())

orderCount <- table(Data$Orders)
orderCount <- orderCount[orderCount > 2]
OrderNames <- names(orderCount)

for (Clade_name in OrderNames) {
  tree <- read.tree("min20Fixed516.nwk")
  filterData <- filter(Data, Orders == Clade_name)
  #for malignancy
  cutData <- filterData[, c(9,17), drop=FALSE]
  # for neo
  #cutData <- filterData[, c(9,13), drop=FALSE]
  cutData$Species <- gsub(" ", "_", cutData$Species) 
  includedSpecies<-cutData$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
  cutData <- cutData[!(cutData$Keep==FALSE),]
  rownames(cutData)<-cutData$Species
  pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10
  x <- setNames(as.numeric(cutData[, 2]), cutData[, 1])
  
  name.check(pruned.tree, x)
  
  results <- tryCatch({
    fitEB_sim <- fitContinuous(pruned.tree, x, model="EB")
    fitBM_sim <- fitContinuous(pruned.tree, x, model="BM")
    fitOU_sim <- fitContinuous(pruned.tree, x, model="OU")
    fitRateTrend_sim <- fitContinuous(pruned.tree, x, model = "rate_trend")
    fitMeanTrend_sim <- fitContinuous(pruned.tree, x, model = "mean_trend")
    fitWhite_sim <- fitContinuous(pruned.tree, x, model = "white")
    
    
    EBresult <- data.frame(Clade = Clade_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x))
    BMresult <- data.frame(Clade = Clade_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x))
    OUresult <- data.frame(Clade = Clade_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x))
    RateTrendresult <- data.frame(Clade = Clade_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x))
    MeanTrendresult <- data.frame(Clade = Clade_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
    Whiteresult<- data.frame(Clade = Clade_name, BestModel = "white", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
    
    testResults <- rbind(EBresult, BMresult, OUresult, RateTrendresult,MeanTrendresult,Whiteresult)
    max_row <- testResults[which.max(testResults$LogLik), ]
    print(max_row)
    resultsClade<-rbind(resultsClade, max_row)
  }, error = function(e) {
    cat("Error encountered: ", e$message, "\nAttempting to reClade tree and retry.\n")
    
    reordered_tree <- reorder.phylo(pruned.tree, order = "postorder")
    
    tryCatch({
      fitEB_sim <- fitContinuous(reCladeed_tree, x, model="EB")
      fitBM_sim <- fitContinuous(reCladeed_tree, x, model="BM")
      fitOU_sim <- fitContinuous(reCladeed_tree, x, model="OU")
      fitRateTrend_sim <- fitContinuous(reCladeed_tree, x, model = "rate_trend")
      fitMeanTrend_sim <- fitContinuous(reCladeed_tree, x, model = "mean_trend")
      fitWhite_sim <- fitContinuous(reCladeed_tree, x, model = "white")
      
      EBresult <- data.frame(Clade = Clade_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x))
      BMresult <- data.frame(Clade = Clade_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x))
      OUresult <- data.frame(Clade = Clade_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x))
      RateTrendresult <- data.frame(Clade = Clade_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x))
      MeanTrendresult <- data.frame(Clade = Clade_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
      Whiteresult<- data.frame(Clade = Clade_name, BestModel = "white", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
      
      
      testResults <- rbind(EBresult, BMresult, OUresult, RateTrendresult,MeanTrendresult,Whiteresult)
      max_row <- testResults[which.max(testResults$LogLik), ]
      resultsClade<- rbind(resultsClade, max_row)
    }, error = function(e) {
      cat("Second attempt failed: ", e$message, "\nReturning NA for all result rows.\n")
      NAresult <- data.frame(Orders = Clade_name, BestModel = NA, LogLik = NA, nSpecies = length(x))
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

ou<-names(speciesModelsOrder)[speciesModelsOrder=="OU"]
lambda<-names(speciesModelsOrder)[speciesModelsOrder=="Lambda"]
EB<-names(speciesModelsOrder)[speciesModelsOrder=="EB"]
#BM<-names(speciesModels)[speciesModels=="BM"]
name.check(pruned.tree, speciesModelsOrder)
reordered_treeClade <- reorder.phylo(pruned.tree, order = "cladewise")




#now for Cladees

Data <- read.csv("min20516.csv")
resultsClade <- data.frame(Clade = character(), BestModel = character(), LogLik = numeric(), nSpecies = numeric())

CladeCount <- table(Data$Clade)
CladeCount <- CladeCount[CladeCount > 2]
CladeNames <- names(CladeCount)

for (Clade_name in CladeNames) {
  tree <- read.tree("min20Fixed516.nwk")
  filterData <- filter(Data, Clade == Clade_name)
  #for malignancy
  cutData <- filterData[, c(9,17), drop=FALSE]
  # for neo
  #cutData <- filterData[, c(9,13), drop=FALSE]
  cutData$Species <- gsub(" ", "_", cutData$Species) 
  includedSpecies<-cutData$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
  cutData <- cutData[!(cutData$Keep==FALSE),]
  rownames(cutData)<-cutData$Species
  pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-10
  x <- setNames(as.numeric(cutData[, 2]), cutData[, 1])
  
  name.check(pruned.tree, x)
  
  results <- tryCatch({
    fitEB_sim <- fitContinuous(pruned.tree, x, model="EB")
    fitBM_sim <- fitContinuous(pruned.tree, x, model="BM")
    fitOU_sim <- fitContinuous(pruned.tree, x, model="OU")
    fitRateTrend_sim <- fitContinuous(pruned.tree, x, model = "rate_trend")
    fitMeanTrend_sim <- fitContinuous(pruned.tree, x, model = "mean_trend")
    fitWhite_sim <- fitContinuous(pruned.tree, x, model = "white")
    
    
    EBresult <- data.frame(Clade = Clade_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x))
    BMresult <- data.frame(Clade = Clade_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x))
    OUresult <- data.frame(Clade = Clade_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x))
    RateTrendresult <- data.frame(Clade = Clade_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x))
    MeanTrendresult <- data.frame(Clade = Clade_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
    Whiteresult<- data.frame(Clade = Clade_name, BestModel = "white", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
    
    testResults <- rbind(EBresult, BMresult, OUresult, RateTrendresult,MeanTrendresult,Whiteresult)
    max_row <- testResults[which.max(testResults$LogLik), ]
    print(max_row)
    resultsClade<-rbind(resultsClade, max_row)
  }, error = function(e) {
    cat("Error encountered: ", e$message, "\nAttempting to reClade tree and retry.\n")
    
    reCladeed_tree <- reorder.phylo(pruned.tree, Clade = "postClade")
    
    tryCatch({
      fitEB_sim <- fitContinuous(reCladeed_tree, x, model="EB")
      fitBM_sim <- fitContinuous(reCladeed_tree, x, model="BM")
      fitOU_sim <- fitContinuous(reCladeed_tree, x, model="OU")
      fitRateTrend_sim <- fitContinuous(reCladeed_tree, x, model = "rate_trend")
      fitMeanTrend_sim <- fitContinuous(reCladeed_tree, x, model = "mean_trend")
      fitWhite_sim <- fitContinuous(reCladeed_tree, x, model = "white")
      
      EBresult <- data.frame(Clade = Clade_name, BestModel = "EB", LogLik = fitEB_sim$opt$lnL, nSpecies = length(x))
      BMresult <- data.frame(Clade = Clade_name, BestModel = "BM", LogLik = fitBM_sim$opt$lnL, nSpecies = length(x))
      OUresult <- data.frame(Clade = Clade_name, BestModel = "OU", LogLik = fitOU_sim$opt$lnL, nSpecies = length(x))
      RateTrendresult <- data.frame(Clade = Clade_name, BestModel = "rate_trend", LogLik = fitRateTrend_sim$opt$lnL, nSpecies = length(x))
      MeanTrendresult <- data.frame(Clade = Clade_name, BestModel = "mean_trend", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
      Whiteresult<- data.frame(Clade = Clade_name, BestModel = "white", LogLik = fitMeanTrend_sim$opt$lnL, nSpecies = length(x))
      
      
      testResults <- rbind(EBresult, BMresult, OUresult, RateTrendresult,MeanTrendresult,Whiteresult)
      max_row <- testResults[which.max(testResults$LogLik), ]
      resultsClade<- rbind(resultsClade, max_row)
    }, error = function(e) {
      cat("Second attempt failed: ", e$message, "\nReturning NA for all result rows.\n")
      NAresult <- data.frame(Clade = Clade_name, BestModel = NA, LogLik = NA, nSpecies = length(x))
      resultsClade<-rbind(resultsClade, NAresult)
    })
  })
}

# Print the final results to check
print(resultsClade)

modelsBestFitClade<-na.omit(resultsClade)






# Filtering the data
filtered_data <- Data %>% 
  filter(Clade %in% as.list(modelsBestFitClade$Clade))

bestFitSpecies<-left_join(filtered_data, modelsBestFitClade, by = "Clade")


cutDataClade<-bestFitSpecies[c(9,5,44)]
tree <- read.tree("min20Fixed516.nwk")

cutDataClade$Species <- gsub(" ", "_", cutDataClade$Species) 
includedSpecies<-cutDataClade$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutDataClade$Keep <- cutDataClade$Species %in% pruned.tree$tip.label
cutDataClade <- cutDataClade[!(cutDataClade$Keep==FALSE),]
rownames(cutDataClade)<-cutDataClade$Species


speciesModels <- setNames((cutDataClade[, 3]), (cutDataClade[, 1]))

ou<-names(speciesModels)[speciesModels=="OU"]
lambda<-names(speciesModels)[speciesModels=="Lambda"]
EB<-names(speciesModels)[speciesModels=="EB"]
#BM<-names(speciesModels)[speciesModels=="BM"]
name.check(pruned.tree, speciesModels)
reordered_tree <- reorder.phylo(pruned.tree, Clade = "cladewise")

cutDataClade <- semi_join(cutDataClade, cutDataOrder, by = "Species")


painted<-reordered_treeClade


# Loop through each row in the dataframe
for (i in 1:nrow(modelsBestFitClade)) {
  clade_data <- filter(cutDataClade, Clade == modelsBestFitClade$Clade[i])
  
  # Check the best model type and apply the corresponding functions
  if (modelsBestFitClade$BestModel[i] == "OU") {
    painted <- paintSubTree(painted, node=findMRCA(painted, clade_data$Species), state="OU", anc="0")
  } else if (modelsBestFitClade$BestModel[i] == "EB") {
    painted <- paintSubTree(painted, node=findMRCA(painted, clade_data$Species), state="EB", anc="0")
  } else if (modelsBestFitClade$BestModel[i] == "Lambda") {
    painted <- paintSubTree(painted, node=findMRCA(painted, clade_data$Species), state="Lambda", anc="0")
  }
  
  
}


# Loop through each row in the dataframe
for (i in 1:nrow(modelsBestFit)) {
  order_data <- filter(cutDataOrder, Orders == modelsBestFit$Orders[i])
  
  # Check the best model type and apply the corresponding functions
  if (modelsBestFit$BestModel[i] == "OU") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="OU", anc="0")
  } else if (modelsBestFit$BestModel[i] == "EB") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="EB", anc="0")
  } else if (modelsBestFit$BestModel[i] == "Lambda") {
    painted <- paintSubTree(painted, node=findMRCA(painted, order_data$Species), state="Lambda", anc="0")
  }
  
  
}





colors<- c("black","#DAA520", "#008080", "#CC5500", "#568203")
names(colors)<-c("0","EB","OU", "rate_trend", "white")
plotSimmap(painted, lwd=2, split.vertical=TRUE, ftype="i", fsize = 0.0000001, colors = colors)



# Loop through each row in the dataframe
for (i in 1:nrow(modelsBestFit)) {
  order_data <- filter(cutDataOrder, Orders == modelsBestFit$Orders[i])
  cladelabels(painted, text =  modelsBestFit$Orders[i], node = findMRCA(painted, order_data$Species), offset = 1, orientation = "horizontal")
  
}






add.simmap.legend(add.simmap.legend(leg=sort(unique(speciesModelsOrder)),
                                    colors=palette()[2:4]))



