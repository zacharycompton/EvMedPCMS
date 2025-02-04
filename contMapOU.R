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


Data <- read.csv("min20516.csv")

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


# fitBM<-fastAnc(tree, x, model = "BM")
# 
# obj <- contMap(tree, x, method = "user", anc.states = fitBM, plot = TRUE)
# obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
# plot(obj, fsize = .8, leg.txt = "Malignancy Prev BM")

# fitBMMCMC_sim<-fitContinuousMCMC(tree, x, model="SSP", Ngens = 10000,
#                                  sampleFreq = 500, printFreq = 10, sample.node.states = TRUE, 
#                                  outputName = "ouMCMC_output")
# 
# OUParams<-read.table("ouMCMC_output_model_params.txt", header = TRUE)
# OUParams$theta
# OUSimNode<-read.table("ouMCMC_output_nodestates.txt", header = TRUE)
# 
# 
# node<-as.data.frame(OUSimNode)
# 
# nodeVal<-tail(node, n=1)
# 
# typeof(pruned.tree$node.label)
# pruned.tree$edge.length
# 
# 
# nodeVals<-c(t(nodeVal))
# cleaned_names <- as.integer(gsub("X|x", "", names(OUSimNode)))
# 
# nodeVals<-setNames(nodeVals,cleaned_names)
# fitBM<-nodeVals
# #fitBM<-fastAnc(tree, x, model= "OU" )

edges <- lambda_tree$edge[lambda_tree$edge[,2] > length(lambda_tree$tip.label), ]

#would like to see if I can get this to work properly with bigger phylo. Check if descendants greater than 2,

branch_info <- data.frame(
  descendant = integer(),  # To store descendant node IDs
  branch_rate = numeric(), # To store calculated branch rates
  stringsAsFactors = FALSE # To avoid converting strings to factors automatically
)

# Loop over edges
for (i in 1:nrow(edges)) {
  
  # Extract ancestor and descendant from the current edge
  ancestor <- edges[i, 1]
  descendant <- edges[i, 2]
  
  # Calculate the absolute change in trait
  trait_change <- abs(as.numeric(pagel[as.character(descendant)]) - as.numeric(pagel[as.character(ancestor)]))
  
  # Check for very small edge lengths and adjust if necessary
  if (lambda_tree$edge.length[i] <= 1e-10) {
    lambda_tree$edge.length[i] = 0
  }
  
  # Get the branch length
  branch_length <- tree$edge.length[i]
  
  # Calculate the branch rate, handling division by zero if needed
  branch_rate <- if (branch_length > 0) trait_change / branch_length else NA
  
  # Append results to the data frame
  branch_info <- rbind(branch_info, data.frame(descendant = descendant, branch_rate = branch_rate))
}
# Associate rates with branches

new_row <- data.frame(descendant = 100, branch_rate = 0, stringsAsFactors = FALSE)

# Prepend the new row to the existing data frame
branch_info <- rbind(new_row, branch_info)


branch_vector <- branch_info$branch_rate

# Set the names of the vector to the 'descendant' column
names(branch_vector) <- branch_info$descendant

# Optional: Ensure the vector is numeric
branch_rates <-branch_vector



branch_rates[is.na(branch_rates)] <- 0

#branch_rates<-na.omit(branch_rates)

branch_rates<-branch_rates[1:lambda_tree$Nnode]

tip_rates <- setNames(vector("numeric", length(lambda_tree$tip.label)), lambda_tree$tip.label)  # Properly named by species

for (tip in 1:length(lambda_tree$tip.label)) {  # Assuming species are nodes 1-98
  # Find the immediate ancestor
  ancestor_index <- which(lambda_tree$edge[,2] == tip)
  if (length(ancestor_index) == 1) {
    ancestor <- lambda_tree$edge[ancestor_index, 1]
    if (!is.na(pagel[as.character(ancestor)]) && !is.na(x[as.character(lambda_tree$tip.label[tip])])) {
      branch_length <- lambda_tree$edge.length[ancestor_index]
      trait_change <- abs(as.numeric(x[as.character(lambda_tree$tip.label[tip])]) - as.numeric(pagel[as.character(ancestor)]))
      tip_rates[as.character(lambda_tree$tip.label[tip])] <- trait_change / branch_length
    } else {
      tip_rates[as.character(lambda_tree$tip.label[tip])] <- NA  # Assign NA if data is missing
    }
  } else {
    tip_rates[as.character(lambda_tree$tip.label[tip])] <- NA  # Assign NA if no single ancestor is found
  }
}



tipsdf<-as.data.frame(tip_rates)
tips<-(tipsdf$tip_rates)
names(tips)<-rownames(tipsdf)
name.check(lambda_tree, tips)

#log is good for seeing HOW much change. No log good for seeing how things change in what direction. 

# # Min-Max Scaling to range [-1, 1]
range_val <- max(abs(tip_rates))
scaled_tip_rates <- tip_rates / range_val

# Min-Max Scaling to range [-1, 1]
range_val <- max(abs(branch_rates))
scaled_branch_rates <- branch_rates / range_val
branch_rates[branch_rates == 0] <- 1e-10
branchdf<-as.data.frame(branch_rates)
branches<-(branchdf$branch_rates)
names(branches)<-rownames(branchdf)




x_binary <- as.integer(x > 0)
names(x_binary)<-names(x)


lambda_tree$edge.length[lambda_tree$edge.length <= 0] <- 1e-10

threshOU<-ancThresh(lambda_tree,x_binary, model = "lambda")


# # Min-Max Scaling to range [-1, 1]
# range_val <- max(abs(branch_rates))
# scaled_branch_rates <- branch_rates / range_val
pallete<-c("darkgrey", "pink")

plot(lambda_tree,no.margin=TRUE,edge.width=2,cex=0.7)
percentage_labels <- sprintf("%.2f%%", branches * 100)
#labelnodes(text = percentage_labels,node=1:tree$Nnode+Ntip(tree),circle.exp = 0.9, interactive = FALSE)
nodelabels(pie=threshOU$ace,piecol=pallete[1:2],cex=0.8)




obj <- contMap(lambda_tree, scaled_tip_rates, method = "user", anc.states = scaled_branch_rates, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj,legend = FALSE, fsize = 0.001, leg.txt = paste("Malignancy Evo Rate Pagel, Lambda ~", round(pagelFit$opt$lambda,2)), type = "fan")
nodelabels(pie=threshOU$ace,piecol=pallete[1:2],cex=0.3)

legend("bottomright",                  # Position of the legend in the plot
       legend = c("No Cancer", "Cancer"),  # Labels for each category
       fill = pallete[1:2],       # Colors used in the pie charts
       title = "Cancer Threshold",      # Title of the legend
       cex = 0.8)    

add.color.bar(100,obj$cols,title= paste("Malignancy Evo Rate Pagel, Lambda ~", round(pagelFit$opt$lambda,2)),
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


#nodelabels(text = percentage_labels,node=1:tree$Nnode+Ntip(tree), bg = "black",cex = 0.8, frame = "circle")

# 
# # Calculate the range of your data focusing on the majority below 2
# small_value_range <- c(min(as.numeric(branches * 100)), 2)
# large_value_range <- c(2, max(as.numeric(branches * 100)))
# 
# # More breaks in the small value range
# small_breaks <- seq(from = small_value_range[1], to = small_value_range[2], length.out = 8)
# 
# # Fewer breaks above 2, to ensure these get differentiated colors
# large_breaks <- seq(from = large_value_range[1], to = large_value_range[2] + 0.1, length.out = 3)  # Adding 0.1 to ensure coverage
# 
# # Combine breaks
# breaks <- c(small_breaks, large_breaks[-1])  # Exclude the first of the large breaks to avoid duplication
# 
# # Define the colors with emphasis on variety below 2
# colors_below_2 <- colorRampPalette(c("blue", "green", "yellow"))(length(small_breaks) - 1)
# colors_above_2 <- c("orange", "red")  # Direct assignment for above 2
# 
# # Combine colors
# colors <- c(colors_below_2, colors_above_2)
# 
# # Assign colors based on intervals
# color_labels <- colors[findInterval(as.numeric(branches * 100), vec = breaks)]
# 
# # Visualization of the tree
# obj <- contMap(tree, x, method = "user", anc.states = fitBM, plot = TRUE)
# obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
# plot(obj, fsize = .8, leg.txt = "Malignancy Prev OU", type = "fan")
# # Applying labels and colors to each node
# for (i in 1:length(percentage_labels)) {
#   nodelabels(text = sprintf("%.2f%%", branches * 100)[i], node = nodeslist[i], cex = 0.3, frame = "circle", bg = color_labels[i])
# }
# add.color.bar(leg=80, cols=colors, title="Rate Change from Parent (%)", lims=c(min(breaks), signif(max(branches * 100)), digits = 2), digits=2, prompt=FALSE, lwd=4, outline=TRUE, x=260, y=-230)
# 
# 
# par("black")
# obj <- contMap(tree, x, method = "user", anc.states = fitBM, plot = TRUE)
# obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
# plot(obj, fsize = .8, leg.txt = "log(Evolutionary Rate) EB")
# labelnodes(text = percentage_labels,node=1:tree$Nnode+Ntip(tree),circle.exp = 0.9, interactive = FALSE)
# 
# 
# 
# 
# fitBMdf<-as.data.frame(fitBM)
# fitBMdf$nodes<-rownames(fitBMdf)
# 
# branchdf$nodes<-rownames(branchdf)
# evoPrevNode<-left_join(branchdf, fitBMdf, by = "nodes")
# colnames(evoPrevNode)<-c("evo", "node", "prev")
# rownames(evoPrevNode)<-evoPrevNode$node
# evoPrevNode$prev<-log(evoPrevNode$prev)
# 
# prevTip<- as.data.frame(x)
# prevTip$tip<-rownames(prevTip)
# 
# Evotips<- as.data.frame(tips)
# Evotips$tip<-rownames(Evotips)
# 
# evoPrevTip<-left_join(prevTip, Evotips, by = "tip")
# colnames(evoPrevTip)<-c("prev", "tip", "evo")
# rownames(evoPrevTip)<-evoPrevTip$tip
# evoPrevTip$prev <- ifelse(evoPrevTip$prev == 0, 0.0001, evoPrevTip$prev)
# evoPrevTip$prev<-log(evoPrevTip$prev)
# 
# 
# evoPrevNode<-evoPrevNode[,c(3,1),drop=FALSE] 
# evoPrevTip<-evoPrevTip[,c(1,3),drop=FALSE] 
# 
# # Standardization (Z-score normalization)
# standardize <- function(x) {
#   return ((x - mean(x)) / sd(x))
# }
# 
# evoPrevNode$evo <- standardize(evoPrevNode$evo)
# evoPrevNode$prev <- standardize(evoPrevNode$prev)
# 
# evoPrevTip$evo <- standardize(evoPrevTip$evo)
# evoPrevTip$prev <- standardize(evoPrevTip$prev)
# 
# # cols<-setNames(rainbow(n =),
# #                levels(as.factor(orderData$Orders)))
# 
# 
# phylomorphospace(tree, evoPrevTip, A=evoPrevNode, label = "off", node.by.map=TRUE,bty="n", cols = cols,las=1 )
# 
