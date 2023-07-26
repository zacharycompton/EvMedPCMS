library(ape)
library(picante)
library(tidyverse)
library(phytools)
library(caper)
library(ggrepel)
library(cowplot)


tree <- read.tree("humanMin20.nwk")
tree
plot(tree)

#functions from pglsSource.R
modPgls.SEy = function (model, data, corClass = corBrownian, tree, se = NULL, 
                        method = c("REML", "ML"), interval = c(0, 1000), corClassValue=1, sig2e=NULL, ...) 
{
  Call <- match.call()
  corfunc <- corClass
  spp <- rownames(data)
  data <- cbind(data, spp)
  if (is.null(se)) 
    se <- setNames(rep(0, Ntip(tree)), tree$tip.label)[spp]
  else se <- se[spp]
  
  lk <- function(sig2e, data, tree, model, ve, corfunc, spp) {
    tree$edge.length <- tree$edge.length * sig2e
    ii <- sapply(1:Ntip(tree), function(x, e) which(e == 
                                                      x), e = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + ve[tree$tip.label]
    vf <- diag(vcv(tree))[spp]
    w <- varFixed(~vf)
    COR <- corfunc(corClassValue, tree, form = ~spp, ...)
    fit <- gls(model, data = cbind(data, vf), correlation = COR, 
               method = method, weights = w)
    -logLik(fit)
  }
  
  if (is.null(sig2e)) {
    fit <- optimize(lk, interval = interval, data = data, tree = tree, 
                    model = model, ve = se^2, corfunc = corfunc, spp = spp)
    sig2e=fit$minimum
  }
  
  tree$edge.length <- tree$edge.length * sig2e
  ii <- sapply(1:Ntip(tree), function(x, e) which(e == x), 
               e = tree$edge[, 2])
  tree$edge.length[ii] <- tree$edge.length[ii] + se[tree$tip.label]^2
  vf <- diag(vcv(tree))[spp]
  w <- varFixed(~vf)
  obj <- gls(model, data = cbind(data, vf), correlation = corfunc(corClassValue, 
                                                                  tree, form = ~spp, ...), weights = w, method = method)
  obj$call <- Call
  obj$sig2e <- sig2e
  obj
}

#Internal function
pglsSEyPagelToOptimizeLambda=function(lambda,model,data,tree,...) {
  -logLik(modPgls.SEy(model=model,data=data,tree=tree,corClassValue=lambda,corClass=corPagel,fixed=T,...)) #Returns -logLikelihood of the pgls.SEy model with lambda fixed to the value of the lambda argument. sig2e will be optimized within modPgls.SEy unless given as an argument here
}

#Function intended for users
pglsSEyPagel=function(model, data, tree, lambdaInterval=c(0,1),...){
  optimizedModel=optimize(pglsSEyPagelToOptimizeLambda,lambdaInterval,model=model,data=data,tree=tree,...) #Optimizes lambda in the lambdaInterval using the pglsSEyPagelToOptimizeLambda function
  return(modPgls.SEy(model=model,data=data,tree=tree,corClass=corPagel,fixed=T,corClassValue=optimizedModel$minimum,...)) #Returns the final model fit
}


## plot tree
plotTree(tree,type="fan", ftype = "i")

## read csv file

primates<-read.csv("min20516.csv",header = TRUE)
#view(primates)

#replaces all space characters in species value of each row to an "_"
primates <- primates %>%
  mutate(Species = gsub(" ", "_", Species))
plot(tree)

## prune the tree to match the data
includedSpecies <- primates$Species
newtips<-str_remove_all(tree$tip.label,"_ott")
newtips<-str_remove_all(newtips,".ott")
newtips<-str_remove_all(newtips,"-ott")
newtips<-str_remove_all(newtips,"[1234567890]")
newtips<-sub('^([^_]+_[^_]+).*', '\\1', newtips)

## pruning the tree
pruned.tree <- tree
tree$tip.label <- newtips
drop.tip(
  pruned.tree, setdiff(
    tree$tip.label, includedSpecies))
keep.tip(pruned.tree,pruned.tree$tip.label)

relate<-cophenetic.phylo(tree)
relate<-as.data.frame(relate)
homoRelate<-cbind(rownames(relate),relate$Homo_sapiens)
colnames(homoRelate) <- c("Species", "relatedness")
view(homoRelate)

primates$Species

cancerRelate<-left_join(primates, homoRelate, by = "Species", copy = TRUE)

#merged_df <- merge(primates, homoRelate, by = "Species", all = TRUE)
#view(merged_df)



#Neo

cutCancer<-cancerRelate[,c(9,10,13,42,43,44),drop=FALSE] 
rownames(cutCancer)<-cutCancer$Species
cutCancer<-na.omit(cutCancer)


cutCancer$relatedness <- as.numeric(cutCancer$relatedness)

view(cutCancer)

SE<-setNames(cutCancer$SE_simple,cutCancer$Species)[rownames(cutCancer)]

relate<-cutCancer$relatedness



pgls.model <- pglsSEyPagel(NeoplasiaPrevalence~relatedness,data=cutCancer,
                           tree=pruned.tree,method="ML",se=SE)
summary(pgls.model)

variable_relate<- data.frame(relatedness = c(0))
hs_pred <- predict(pgls.model, newdata = variable_relate)

hs_pred[1]


#Mal


cutMal<-cancerRelate[,c(9,10,17,42,43,44),drop=FALSE] 
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





