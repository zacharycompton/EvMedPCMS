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

pagel<-pagel$ace

# 
#  pagel<-vector(); 
#  N<-length(lambda_tree$tip); M<-lambda_tree$Nnode
#  
#   
#  # compute the PIC "root" state for each internal node
#    for(i in 1:M+N){
#      pagel[i-N]<-ace(x,multi2di(root(pruned.tree,node=i)),
#                  method="pic")$ace[1]
#     names(pagel)[i-N]<-i
#   }


obj <- contMap(lambda_tree, x, method = "user", anc.states = pagel, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj, fsize = .8, leg.txt = "Malignancy Prev Pagel")




x_binary <- as.integer(x > 0)
names(x_binary)<-names(x)


lambda_tree$edge.length[lambda_tree$edge.length <= 0] <- 1e-10

#threshOU<-ancThresh(lambda_tree,x_binary, model = "lambda")


# # Min-Max Scaling to range [-1, 1]
# range_val <- max(abs(branch_rates))
# scaled_branch_rates <- branch_rates / range_val
pallete<-c("black", "white")

plot(lambda_tree,no.margin=TRUE,edge.width=2,cex=0.7)
#percentage_labels <- sprintf("%.2f%%", branches * 100)
#labelnodes(text = percentage_labels,node=1:tree$Nnode+Ntip(tree),circle.exp = 0.9, interactive = FALSE)
nodelabels(pie=pagel,cex=0.8)




obj <- contMap(lambda_tree, x, method = "user", anc.states = pagel, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj,legend = FALSE, fsize = 0.001, leg.txt = paste("Malignancy Evo Rate Pagel, Lambda ~", round(pagelFit$opt$lambda,2)), type = "fan")
nodelabels(pie=pagel,piecol=pallete[1:2],cex=0.3)

legend("bottomright",                  # Position of the legend in the plot
       legend = c("Malignancy Prevalence"),  # Labels for each category
       fill = pallete[1:1],       # Colors used in the pie charts
       cex = 0.8)    

add.color.bar(100,obj$cols,title= paste("Malignancy Prevalence, Lambda ~", round(pagelFit$opt$lambda,2)),
              lims=obj$lims,digits=3,prompt=FALSE,x=-150,
              y=-200,lwd=2,fsize=0.8,subtitle="")



rodentia<-filter(cutDataAll, is.element(Orders, c("Rodentia")))
primates<-filter(cutDataAll, is.element(Orders, c("Primates")))
cetacea<-filter(cutDataAll, is.element(Orders, c("Cetacea")))
artio<-filter(cutDataAll, is.element(Orders, c("Artiodactyla")))

cetartio<-rbind(cetacea,artio)
carnivora<-filter(cutDataAll, is.element(Orders, c("Carnivora")))
diprotodontia<-filter(cutDataAll, is.element(Orders, c("Diprotodontia")))
chiroptera<-filter(cutDataAll, is.element(Orders, c("Chiroptera")))

aves<-filter(Data, Class == "Aves")
par(fg="#222222")
arc.cladelabels(text="Rodentia",cex = 1.05,node=findMRCA(lambda_tree, rodentia$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)

par(fg="#222222")
arc.cladelabels(text="Cetartiodactyla",cex = 1.05,node=findMRCA(lambda_tree, cetartio$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)
par(fg="#222222")
arc.cladelabels(text="Carnivora",cex = 1.05,node=findMRCA(lambda_tree, carnivora$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)


par(fg="#222222")
arc.cladelabels(text=" Chiroptera",cex = 1.05,node=findMRCA(lambda_tree, chiroptera$Species),ln.offset=1.05,lab.offset=1.1,mark.node =FALSE)

par(fg="#222222")
arc.cladelabels(text="Primates",cex = 1.05,node=findMRCA(lambda_tree, primates$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)

par(fg="#222222")
arc.cladelabels(text="Diprotodontia",cex = 1.05,node=findMRCA(lambda_tree, diprotodontia$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)

