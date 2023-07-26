library(phytools)
library(tidyverse)
packageVersion("phytools")

#read data
Data <- read.csv("discretemin20516.csv")

Data<- filter(Data, is.element(Clade, c("Mammalia")))

tree <- read.tree("min20Fixed516.nwk")
tree<-untangle(ladderize(tree),"read.tree")


#cut data

cutData <- Data
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



cutData$discreteVar <- ifelse(cutData$discreteVar>= 4, 3,cutData$discreteVar)


#standard Rate
#this is where you need to start on your own
#after the $ you need to insert the last column you created

cancerRisk<-setNames(cutData$discreteVar,rownames(cutData))

levels(cancerRisk)<- c("0","<=5%","<=10%",">10%")
head(cancerRisk)


ordered_Mk<-fitHRM(pruned.tree,cancerRisk,niter=10,
                   parallel=TRUE,ncores=10,umbral=TRUE,ordered=TRUE,
                   order=sort(unique(cancerRisk)),pi="fitzjohn",ncat=1,
                   logscale=TRUE)



ordered_Mk

plot(as.Qmatrix(ordered_Mk),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)


plot(ordered_Mk,spacer=0.3,offset=0.03,
     mar=rep(0.1,4))


risk_ancr<-ancr(ordered_Mk,tips=TRUE)
risk_ancr




h<-max(nodeHeights(pruned.tree))
plotTree(pruned.tree,ftype="off",lwd=1,
         ylim=c(0,1.05*h),direction="upwards")
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(c("#f0ead6","darkred","lightblue", "lavender"),
               colnames(risk_ancr$ace))
for(i in 1:Ntip(pruned.tree)){
  tcp<-c(0,cumsum(risk_ancr$ace[pruned.tree$tip.label[i],]))
  tcp<-tcp*0.05*h
  for(j in 2:length(tcp)){
    polygon(rep(pp$xx[i],4)+c(-0.5,-0.5,0.5,0.5),
            h+c(tcp[j-1],tcp[j],tcp[j],tcp[j-1]),
            border=FALSE,col=cols[j-1])
  }
}
legend(x=0,y=0.5*h,c("0","<=5%","<=10%",">10%"),pch=15,
       col=cols,pt.cex=1.5,cex=0.8,bty="n")
par(fg="transparent")
nodelabels(pie=risk_ancr$ace[1:Nnode(pruned.tree)+
                                 Ntip(pruned.tree),],piecol=cols,cex=0.3)




#hidden rates


hrm.umbral<-fitHRM(pruned.tree,cancerRisk,umbral=TRUE,
                   pi="fitzjohn",opt.method="optimParallel",rand_start=TRUE)



hrm.umbral


plot(hrm.umbral,spacer=0.3,offset=0.03,
     mar=rep(0.1,4))

plot(as.Qmatrix(hrm.umbral),width=TRUE,color=TRUE,xlim=c(-1.5,1),
     text=FALSE,max.lwd=4)


risk_ancr<-ancr(hrm.umbral,tips=TRUE)
risk_ancr




h<-max(nodeHeights(pruned.tree))
plotTree(pruned.tree,ftype="off",lwd=1,
         ylim=c(0,1.05*h),direction="upwards")
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(c("#f0ead6","#fffff2","darkred","blue","lightblue","purple","lavender", "#FFD580"),
               colnames(risk_ancr$ace))
for(i in 1:Ntip(pruned.tree)){
  tcp<-c(0,cumsum(risk_ancr$ace[pruned.tree$tip.label[i],]))
  tcp<-tcp*0.05*h
  for(j in 2:length(tcp)){
    polygon(rep(pp$xx[i],4)+c(-0.5,-0.5,0.5,0.5),
            h+c(tcp[j-1],tcp[j],tcp[j],tcp[j-1]),
            border=FALSE,col=cols[j-1])
  }
}
legend(x=0,y=0.5*h,c("0","0*","<=5%","<=5%*","<=10%","<=10%*",">10%",">10%*"),pch=15,
       col=cols,pt.cex=1.5,cex=0.8,bty="n")
par(fg="transparent")
nodelabels(pie=risk_ancr$ace[1:Nnode(pruned.tree)+
                               Ntip(pruned.tree),],piecol=cols,cex=0.3)





nodeVal<-risk_ancr$ace
pruned.tree
##
nodeVal_1 <- nodeVal[-c(1:87),-c(1,3,5,7)]
tipVal<-nodeVal[c(1:87),-c(1,3,5,7)]
hiddensum <- c(rowSums(nodeVal_1, na.rm = FALSE, dims = 1L))
hiddensumTip <- c(rowSums(tipVal, na.rm = FALSE, dims = 1L))


typeof(hiddensum)

#^this should say "double


#make a list of species names and NeoplasiaPrevalence and assign it to variable name "x"
x<-hiddensumTip
#now we have a tree
par(fg="black")
obj<-contMap(pruned.tree,x,plot=TRUE,method="user",anc.states=hiddensum, fsize = .5 )


## invert color scale so red is highest trait value and blue is lowest
obj<-setMap(obj,invert=TRUE)
plot(obj, fsize= .5)

obj

