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

## simulate Brownian evolution on a tree with fastBM
x<-fastBM(pruned.tree,internal=TRUE)

x_ou<-fastBM(pruned.tree,a=0,theta=3,alpha=0.2,sig2=0.1, internal=TRUE)



bmDF<-data.frame(x=x[1:229], mass=cutData$adult_weight.g.)

brownLm<-lm(x~mass,bmDF )

summary(brownLm)




ouDF<-data.frame(x= x_ou[1:229], mass = cutData$adult_weight.g.)


ouLm<- lm(x~mass, ouDF)

summary(ouLm)



martin_cor <- sim.corrs(pruned.tree, corMartins(1, pruned.tree))


