library(phytools)
library(geiger)
library(tidyverse)
library(cowplot)
library(mvMORPH)
library(RevGadgets)
#install.packages("RColorBrewer")
library(ratematrix)
library(RRphylo)
#library(RColorBrewer)
## simulate data under model="BM"

#x<-fastBM(tree)
#plotTree(tree)



#read data
Data <- read.csv("min20516.csv")
Data<- filter(Data, is.element(Orders, c("Primates")))
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
ouData<-x
rownames(ouData)<-ouData$Species

all<-c(t(x))
species<-c(t(x$Species))
neo<-as.numeric(c(t(x$NeoplasiaPrevalence)))

tree<-pruned.tree

x <- setNames(neo, species)

length(x)

typeof(x)

name.check(pruned.tree, ouData)


fitEB_sim <- fitContinuous(tree, x, model="EB")
fitBM_sim <- fitContinuous(tree, x, model="BM")
fitOU_sim <- fitContinuous(tree, x, model="OU")
fitLambda_sim <- fitContinuous(tree, x, model="lambda")
fitLambda_sim <- fitContinuous(tree, x, model="lambda")

fitEB_sim$

r<-RRphylo(tree,x)
r$rates
r$aces
nodeRate<-head(r$rates, n= 97)
nodeRates<-setNames(nodeRate[,1],seq(99,195))
typeof(nodeRates)

tipRate<-tail(r$rates, n=98 )
tipRates<-setNames(tipRate[,1], species)


tree<-pruned.tree

x <- setNames(neo, species)


name.check(tree, tipRate)

length(tree$tip.label)
length(tipRate)
length(na.omit(tipRate))


# Define the endpoints for the color palette (dark blue, green, red)
colors <- c("#00008B","#00a5ff","#00e1ad" ,"#FFFF00", "#FFA500","#FF0000")

# Create the color palette function with 20 colors
color_palette <- colorRampPalette(colors)(2000)




obj<-contMap(tree,tipRates,plot=TRUE,method="user",anc.states=nodeRates)

plot(setMap(obj,invert=TRUE), fsize = .8, leg.txt = "Evolutionary Rate")
mammal.contMap<-setMap(obj,colors)
                                

plot(mammal.contMap,fsize=c(0.7,0.8),
     leg.txt="Evo Rate")


#OU reconstruction

OUrecon<-anc.ML(tree, tipRates, model = "OU")

objOU<-contMap(tree,tipRates,plot=TRUE,method="user",anc.states=OUrecon$ace)

mammal.contMapOU<-setMap(objOU,color_palette)


plot(mammal.contMapOU,fsize=c(0.7,0.8),
     leg.txt="Evo Rate")

ouData<-ouData[,c(2),drop=FALSE] 

OUwie(pruned.tree, ouData, model= "OU1")

obj

par(fg="white")
plotSimmap(obj$tree,obj$cols,type="phylogram", fsize =.00001 )

mammal.contMap<-setMap(mammal.contMap,
                       c("white","#FFFFB2","#FECC5C","#FD8D3C",
                                "#E31A1C"))
add.color.bar(100,cols=obj$cols,title="BM Neoplasia Prevalence", lims=range(x),digits=2)


