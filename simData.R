library(phytools)
library(geiger)
library(tidyverse)
library(cowplot)
library(mvMORPH)
library(RevGadgets)

## simulate data under model="BM"

#x<-fastBM(tree)
#plotTree(tree)


  
  #read data
  Data <- read.csv("min20516.csv")
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
  
  typeof(x)
  
  name.check(tree, x)


fitEB_sim <- fitContinuous(tree, x, model="EB")
fitBM_sim <- fitContinuous(tree, x, model="BM")
fitOU_sim <- fitContinuous(tree, x, model="OU")

evol_rate <- ace(x, tree)


aicScores[nrow(aicScores) + 1,] = c(AIC(fitOU_sim), AIC(fitBM_sim),AIC(fitEB_sim))


aicScores




fitOUMCMC_sim<-fitContinuousMCMC(tree, x, model="SSP", Ngens = 10000,
                  sampleFreq = 500, printFreq = 10, sample.node.states = TRUE, 
                  outputName = "ouMCMC_output")
summary(fitOUMCMC_sim)
fitBMMCMC_sim<-fitContinuousMCMC(tree, x, model="BM", Ngens = 10000,
                  sampleFreq = 500, printFreq = 10, sample.node.states = TRUE,
                  outputName = "bmMCMC_output")

fitEBMCMC_sim<-fitContinuousMCMC(tree, x, model="ACDC.exp", Ngens = 10000,
                  sampleFreq = 500, printFreq = 10,sample.node.states = FALSE,
                  outputName = "ebMCMC_output")

OUsimtxt<-read.table("ouMCMC_output_model_params.txt", header = TRUE)


ggplot(data = OUsimtxt,aes(x= generation,y=logLk ))+
  geom_jitter()+
  scale_y_continuous(
    breaks = c(round(min(OUsimtxt$logLk)),round(median(OUsimtxt$logLk)), round(max(OUsimtxt$logLk)) ),
    labels = c(round(min(OUsimtxt$logLk)),round(median(OUsimtxt$logLk)), round(max(OUsimtxt$logLk))))+
  theme_cowplot(12)+
  geom_line()+
  labs(title = "Log Likelihood per Generation of MCMC OU")


BMsimtxt<-read.table("bmMCMC_output_model_params.txt", header = TRUE)


ggplot(data = BMsimtxt,aes(x= generation,y=logLk ))+
  geom_jitter()+
  scale_y_continuous(
    breaks = c(round(min(BMsimtxt$logLk)),round(median(BMsimtxt$logLk)), round(max(BMsimtxt$logLk)) ),
    labels = c(round(min(BMsimtxt$logLk)),round(median(BMsimtxt$logLk)), round(max(BMsimtxt$logLk))))+
  theme_cowplot(12)+
  geom_line()+
  labs(title = "Log Likelihood per Generation of MCMC BM")


EBsimtxt<-read.table("ebMCMC_output_model_params.txt", header = TRUE)


ggplot(data = EBsimtxt,aes(x= generation,y=logLk ))+
  geom_jitter()+
  theme_cowplot(12)+
  geom_line()+
  labs(title = "Log Likelihood per Generation of MCMC EB")

##this is a frankenstein of heatmap tree and labeling function should work just fine




#subset

specData <- Data[,c(5,9,13),drop=FALSE] 
specData <- na.omit(specData)
specData$Species <- gsub(" ", "_", specData$Species) 
includedSpecies<-specData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
specData$Keep <- specData$Species %in% pruned.tree$tip.label
specData <- specData[!(specData$Keep==FALSE),]
pruned.tree$edge.length[which(pruned.tree$edge.length <=0)] <- 1e-10



mammals<-filter(specData, is.element(Clade, c("Mammalia")))
saur<-filter(specData, is.element(Clade, c("Sauropsida")))
amph<-filter(specData, is.element(Clade, c("Amphibia")))

#plot with labels



OUSimNode<-read.table("ouMCMC_output_nodestates.txt", header = TRUE)

node<-as.data.frame(OUSimNode)

nodeVal<-tail(node, n=1)

typeof(pruned.tree$node.label)
pruned.tree$edge.length


nodeVals<-c(t(nodeVal))
nodeVals<-setNames(nodeVals,seq(293,583))

typeof(nodeVals)



obj<-contMap(pruned.tree,x,plot=TRUE,method="user",anc.states=nodeVals)


## invert color scale so red is highest trait value and blue is lowest
obj<-setMap(obj,invert=TRUE)


obj
par(fg="white")

plotSimmap(obj$tree,obj$cols,type="fan", fsize = .00001 )
add.color.bar(100,cols=obj$cols,title="OU Neoplasia Prevalence", lims=range(x),digits=2)

par(fg="#3B4992FF")
arc.cladelabels(text="Mammalia",cex = .8,node=findMRCA(pruned.tree, mammals$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#EE0000FF")
arc.cladelabels(text="Sauropsida",cex = .8,node=findMRCA(pruned.tree, saur$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#008280FF")
arc.cladelabels(text="Amphibia",cex = .8,node=findMRCA(pruned.tree, amph$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)






#brownian
bmSimNode<-read.table("bmMCMC_output_nodestates.txt", header = TRUE)

node<-as.data.frame(bmSimNode)

nodeVal<-tail(node, n=1)

typeof(pruned.tree$node.label)

nodeVals<-c(t(nodeVal))
nodeVals<-setNames(nodeVals,seq(293,583))


typeof(nodeVals)



obj<-contMap(pruned.tree,x,plot=TRUE,method="user",anc.states=nodeVals)


## invert color scale so red is highest trait value and blue is lowest
obj<-setMap(obj,invert=TRUE)


obj

par(fg="white")
plotSimmap(obj$tree,obj$cols,type="fan", fsize =.00001 )
add.color.bar(100,cols=obj$cols,title="BM Neoplasia Prevalence", lims=range(x),digits=2)

par(fg="#3B4992FF")
arc.cladelabels(text="Mammalia",cex = .8,node=findMRCA(pruned.tree, mammals$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#EE0000FF")
arc.cladelabels(text="Sauropsida",cex = .8,node=findMRCA(pruned.tree, saur$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)

par(fg="#008280FF")
arc.cladelabels(text="Amphibia",cex = .8,node=findMRCA(pruned.tree, amph$Species),ln.offset=1.05,lab.offset=1.15,mark.node=FALSE)





