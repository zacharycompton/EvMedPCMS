library(phytools)
library(geiger)
library(tidyverse)
library(cowplot)
library(mvMORPH)
library(RevGadgets)
#install.packages("RColorBrewer")
library(ratematrix)
library(OUwie)
library(RRphylo)
#library(RColorBrewer)


Data <- read.csv("min20-2022.05.16.csv")

tree <- read.tree("min20Fixed516.nwk")
#for malignancy

Data<-filter(Data, Class == "Mammalia")

#Data<-filter(Data, Orders == "Primates")


cutDataAll <- Data[, c(9,17, 6), drop=FALSE]
cutDataAll$Species <- gsub(" ", "_", cutDataAll$Species) 




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


tree<-multi2di(tree,random=TRUE)
is.binary(tree)

name.check(pruned.tree, x)

#edgelabels()
N<-length(pruned.tree$tip.label)

#fit<-fitContinuous(tree, x, model = "OU")

pagelFit<- fitContinuous(pruned.tree, x, model="lambda")

lambda_tree <- rescale(pruned.tree, model = "lambda", lambda = pagelFit$opt$lambda)
lambda_tree$edge.length[lambda_tree$edge.length <= 0] <- 1e-10

name.check(lambda_tree, x)



pagel<-anc.ML(lambda_tree, x, model = "BM")

A<-fastAnc(lambda_tree,x,CI=TRUE)





tree<-paintSubTree(pruned.tree,node=length(pruned.tree$tip)+1,"1")
# our transparencies
trans <- as.character(floor(0:50 / 2))
trans[as.numeric(trans) < 10] <- paste("0", trans[as.numeric(trans) < 10], sep="")

# plot
for(i in 0:50){
  p<-i/length(trans)
  phenogram(tree,c(x,(1-p)*A$CI95[,1]+p*A$ace), colors=setNames(paste("#0000ff",trans[i+1],sep=""),1), add=i>0, fsize=0, ylab = "Malignancy Prevalence", xlab = "Time (MYA)", spread.labels = FALSE)
  phenogram(tree,c(x,(1-p)*A$CI95[,2]+p*A$ace), colors=setNames(paste("#0000ff",trans[i+1],sep=""),1), add=TRUE, fsize=0, ylab = "Malignancy Prevalence", xlab = "Time (MYA)", spread.labels = FALSE)
}
phenogram(tree, c(x, A$ace), add = TRUE, colors = setNames("white", 1), ylab = "Malignancy Prevalence", xlab = "Time (MYA)", fsize = 0, spread.labels = FALSE)

