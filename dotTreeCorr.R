library(phytools)
library(geiger)
library(tidyverse)
library(cowplot)
library(mvMORPH)
library(RevGadgets)
library(RRphylo)

remotes::install_github("Rekyt/divr")


library(divr)




#read data
Data <- read.csv("min20516.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
View(Data)


tree <- read.tree("min20Fixed516.nwk")


#cut data

cutData <- Data[,c(9,13,17, 30,38,40),drop=FALSE] 
cutData$Species <- gsub("NaN", "NA", cutData$Species) 
cutData[cutData == -1 ] <-NA
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

rownames(cutData)<-cutData$Species
cutData<-cutData[,c(2:6),drop=FALSE] 

name.check(pruned.tree, cutData)


#For Matt: Pearson Correlation


data <- data.frame(lapply(cutData, as.numeric))

rownames(data)<-rownames(cutData)

mtxData<-as.matrix(cutData)


vcv<-vcv.phylo(pruned.tree)

phyCorr<-phyl.vcv(mtxData,vcv(pruned.tree),1)

#Matt: Refer to this data object to make dot tree

phyCorr$R

#transform

transData<-cutData
transData$adult_weight.g.<-log(transData$adult_weight.g.)
transData$max_longevity.months.<-log(transData$max_longevity.months.)

transData$NeoplasiaPrevalence<-transData$NeoplasiaPrevalence*100
transData$MalignancyPrevalence<-transData$MalignancyPrevalence*100

name.check(pruned.tree, transData)






dotTree(pruned.tree,transData,fsize=0.5,border="transparent",cex.dot=1.5,
        length=20,standardize=TRUE,colors=palette()[2], legend = FALSE)


palette()[2]

