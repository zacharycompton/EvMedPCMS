library(phytools)
library(ape)
library(dplyr)
library(OUwie)




#read data
Data <- read.csv("min20516.csv")
View(Data)





#adult weight models
#adult weight neo

cutData <- Data[,c(5,9,10,11,13,38,42),drop=FALSE] 
cutData[cutData$adult_weight == -1, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]


# Simulate Brownian Motion values on the tree for two variables
bm_sim_x <- fastBM(pruned.tree, mu= 1)
bm_sim_y <- fastBM(pruned.tree, mu = 1)



bmLinear<-lm(bm_sim_y~bm_sim_x)

summary(bmLinear)


# Simulate OU values on the tree for two variables
ou_sim_x<-fastBM(pruned.tree,a=0,theta=3,alpha=0.2,sig2=0.1)
ou_sim_y<-fastBM(pruned.tree,a=0,theta=3,alpha=0.2,sig2=0.1)



ouLinear<-lm(ou_sim_y~ou_sim_x)

summary(ouLinear)

simMartins <- function(tree, alpha, sig2, theta, a0, internal) {
  tree <- reorder(tree, "cladewise")
  num_traits <- 1  # Adjust this to the desired number of traits
  X <- matrix(0, nrow = length(tree$edge), ncol = num_traits + 1)
  root <- length(tree$tip.label) + 1
  X[which(tree$edge[, 1] == root), 1] <- a0
  
  for (trait in 2:(num_traits + 1)) {
    for (i in 1:(nrow(X) - 1)) {
      t <- tree$edge.length[i]
      s2 <- sig2 * (1 - exp(-2 * alpha * t)) / (2 * alpha)
      X[i + 1, trait] <- exp(-alpha * t) * X[i, trait] + (1 - exp(-alpha * t)) * theta + rnorm(n = 1, sd = sqrt(s2))
    }
  }
  
  x <- sapply(1:max(tree$edge), function(x, y, tree) y[tree$edge[, 2] == x, trait], y = X, tree = tree)
  x <- setNames(x, c(tree$tip.label, 1:tree$Nnode + length(tree$tip.label)))
  
  if (internal == TRUE) {
    return(x)  # include internal nodes
  } else {
    return(x[1:length(tree$tip.label)])  # tip nodes only
  }
}



# Simulate OU values on the tree for two variables
marty_sim_y<-simMartins(pruned.tree,a=0,theta=3,alpha=0.2,sig2=0.1, internal = FALSE)
marty_sim_x<-simMartins(pruned.tree,a=0,theta=3,alpha=0.2,sig2=0.1, internal = FALSE)


valuesSimY <- sapply(marty_sim_y, `[`, 1)
valuesSimX <- sapply(marty_sim_x, `[`, 1)


martyLinear<-lm(valuesSimY~valuesSimX)

summary(martyLinear)










