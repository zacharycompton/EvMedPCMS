library(phytools)
library(geiger)
library(tidyverse)
library(cowplot)
library(mvMORPH)
library(RevGadgets)
#install.packages("RRphylo")
library(RRphylo)
## simulate data under model="BM"

#x<-fastBM(tree)
#plotTree(tree)



#read data
Data <- read.csv("min20516.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
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

length(x)

typeof(x)

name.check(tree, x)


fitEB_sim <- fitContinuous(tree, x, model="EB")
fitBM_sim <- fitContinuous(tree, x, model="BM")
fitOU_sim <- fitContinuous(tree, x, model="OU")

r<-RRphylo(tree,x)
r$rates
nodeRate<-head(r$rates, n= 97)
nodeRates<-setNames(nodeRate[,1],seq(99,195))
typeof(nodeRates)

tipRate<-tail(r$rates, n=98 )
length(tipRate)

tipRates<-setNames(tipRate[,1], species)



obj<-contMap(tree,tipRates,plot=TRUE,method="user",anc.states=nodeRates)

plot(setMap(obj,invert=TRUE))
mammal.contMap<-setMap(obj,
                       c("#0000ff", "#0078ff", "#00a5ff", "#00c6d9", "#00e1ad", "#46e085", "#74dd54", "#a0d600" ,
                                  "#c2b700", "#de9200", "#f36300", "#ff0000"))
                                

plot(mammal.contMap,fsize=c(0.7,0.8),
     leg.txt="Evo Rate")

obj

par(fg="white")
plotSimmap(obj$tree,obj$cols,type="phylogram", fsize =.00001 )

mammal.contMap<-setMap(mammal.contMap,
                       c("white","#FFFFB2","#FECC5C","#FD8D3C",
                                "#E31A1C"))
add.color.bar(100,cols=obj$cols,title="BM Neoplasia Prevalence", lims=range(x),digits=2)


