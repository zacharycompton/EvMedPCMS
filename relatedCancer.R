library(ape)
library(picante)
library(tidyverse)
library(phytools)
library(caper)
library(ggrepel)
library(cowplot)


tree <- read.tree("homoMonky.nwk")
tree
plot(tree)

## plot tree
plotTree(tree,type="fan", ftype = "i")

## read csv file

primates<-read.csv("min20Primates.csv",header = TRUE)  
row.names(primates)= gsub(" ", "_", row.names(primates))

plot(tree)

## prune the tree to match the data
includedSpecies <- primates$Species
newtips<-str_remove_all(tree$tip.label,"_ott")
newtips<-str_remove_all(newtips,".ott")
newtips<-str_remove_all(newtips,"-ott")
newtips<-str_remove_all(newtips,"[1234567890]")
newtips<-sub('^([^_]+_[^_]+).*', '\\1', newtips)

## pruning the tree
tree$tip.label <- newtips
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)

relate<-cophenetic.phylo(tree)
relate<-as.data.frame(relate)
homoRelate<-cbind(rownames(relate),relate$Homo_sapiens)
view(homoRelate)
colnames(homoRelate) <- c("Species", "relatedness")

primates$Species

cancerRelate<-left_join(primates, homoRelate, by = "Species", copy = TRUE)

merged_df <- merge(primates, homoRelate, by = "species", all = TRUE)

view(cancerRelate)


#Neo

cutCancer<-cancerRelate[,c(8,9,12,46,48),drop=FALSE] 
rownames(cutCancer)<-cutCancer$Species
cutCancer<-na.omit(cutCancer)

cutCancer$relatedness <- as.numeric(cutCancer$relatedness)

SE<-setNames(cutCancer$SE_simple,cutCancer$Species)[rownames(cutCancer)]

relate<-cutCancer$relatedness



pgls.model <- pglsSEyPagel(NeoplasiaPrevalence~relatedness,data=cutCancer,
                           tree=pruned.tree,method="ML",se=SE)
summary(pgls.model)

variable_relate<- data.frame(relatedness = c(0))
hs_pred <- predict(pgls.model, newdata = variable_relate)

hs_pred[1]


#Mal


cutMal<-cancerRelate[,c(8,9,15,46,48),drop=FALSE] 
rownames(cutMal)<-cutMal$Species
cutMal<-na.omit(cutMal)

cutMal$relatedness <- as.numeric(cutMal$relatedness)

SE<-setNames(cutMal$SE_simple,cutMal$Species)[rownames(cutMal)]

relate<-cutMal$relatedness



pgls.model <- pglsSEyPagel(MalignancyPrevalence~relatedness,data=cutMal,
                           tree=pruned.tree,method="ML",se=SE)
summary(pgls.model)

variable_relate<- data.frame(relatedness = c(0))
hs_pred <- predict(pgls.model, newdata = variable_relate)

hs_pred[1]





