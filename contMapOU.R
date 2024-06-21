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


#read data
Data <- read.csv("min20516.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
View(Data)


tree <- read.tree("min20Fixed516.nwk")


#cut data

cutData <- Data[,c(9,13,6),drop=FALSE] 
cutData <- na.omit(cutData)
cutData$Species <- gsub(" ", "_", cutData$Species) 
#orderData<-cutData
#orderData<-left_join(ouData, orderData, by = "Species")
#orderData<-orderData[,c(1,4),drop=FALSE] 
#rownames(orderData)<-orderData$Species
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
x<-cutData[,c(1,2),drop=FALSE] 
ouData<-x
rownames(ouData)<-cutData$Species

ouData<-ouData[,c(2),drop=FALSE] 
pruned.tree$edge.length[which(pruned.tree$edge.length <=0)] <- 1e-10


data(tworegime)
trait

OUwie(tree,trait,model=c("OUMV"))


all<-c(t(x))
species<-c(t(x$Species))
neo<-as.numeric(c(t(x$NeoplasiaPrevalence)))

tree<-pruned.tree

x <- setNames(neo, species)

tree<-multi2di(tree,random=TRUE)
is.binary(tree)

name.check(tree, x)

#edgelabels()
N<-length(pruned.tree$tip.label)

#fit<-fitContinuous(tree, x, model = "OU")


fitOU<-reconstruct(x,tree, method = "GLS_OUS")
fitOU<-fitOU$ace


obj <- contMap(tree, x, method = "user", anc.states = fitOU, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj, fsize = .8, leg.txt = "Neoplasia Prev OU")


# fitBM<-fastAnc(tree, x, model = "BM")
# 
# obj <- contMap(tree, x, method = "user", anc.states = fitBM, plot = TRUE)
# obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
# plot(obj, fsize = .8, leg.txt = "Neoplasia Prev BM")

# fitOUMCMC_sim<-fitContinuousMCMC(tree, x, model="SSP", Ngens = 10000,
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
# fitOU<-nodeVals
# #fitOU<-fastAnc(tree, x, model= "OU" )

edges <- tree$edge[tree$edge[,2] > length(tree$tip.label), ]

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
  trait_change <- abs(as.numeric(fitOU[as.character(descendant)]) - as.numeric(fitOU[as.character(ancestor)]))
  
  # Check for very small edge lengths and adjust if necessary
  if (tree$edge.length[i] <= 1e-10) {
    tree$edge.length[i] = 0
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

branch_rates<-branch_rates[1:tree$Nnode]

tip_rates <- setNames(vector("numeric", length(tree$tip.label)), tree$tip.label)  # Properly named by species

for (tip in 1:length(tree$tip.label)) {  # Assuming species are nodes 1-98
  # Find the immediate ancestor
  ancestor_index <- which(tree$edge[,2] == tip)
  if (length(ancestor_index) == 1) {
    ancestor <- tree$edge[ancestor_index, 1]
    if (!is.na(fitOU[as.character(ancestor)]) && !is.na(x[as.character(tree$tip.label[tip])])) {
      branch_length <- tree$edge.length[ancestor_index]
      trait_change <- abs(as.numeric(x[as.character(tree$tip.label[tip])]) - as.numeric(fitOU[as.character(ancestor)]))
      tip_rates[as.character(tree$tip.label[tip])] <- trait_change / branch_length
    } else {
      tip_rates[as.character(tree$tip.label[tip])] <- NA  # Assign NA if data is missing
    }
  } else {
    tip_rates[as.character(tree$tip.label[tip])] <- NA  # Assign NA if no single ancestor is found
  }
}



tipsdf<-as.data.frame(tip_rates)
tips<-(tipsdf$tip_rates)
names(tips)<-rownames(tipsdf)
name.check(tree, tips)

#log is good for seeing HOW much change. No log good for seeing how things change in what direction. 

# # Min-Max Scaling to range [-1, 1]
# range_val <- max(abs(tip_rates))
# scaled_tip_rates <- tip_rates / range_val

# Min-Max Scaling to range [-1, 1]
#range_val <- max(abs(branch_rates))
#scaled_branch_rates <- branch_rates / range_val
branch_rates[branch_rates == 0] <- 1e-10
branchdf<-as.data.frame(branch_rates)
branches<-(branchdf$branch_rates)
names(branches)<-rownames(branchdf)

# # Min-Max Scaling to range [-1, 1]
# range_val <- max(abs(branch_rates))
# scaled_branch_rates <- branch_rates / range_val

plot(tree,no.margin=TRUE,edge.width=2,cex=0.7)
percentage_labels <- sprintf("%.2f%%", branches * 100)
labelnodes(text = percentage_labels,node=1:tree$Nnode+Ntip(tree),circle.exp = 0.9, interactive = FALSE)


obj <- contMap(tree, x, method = "user", anc.states = fitOU, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj, fsize = .8, leg.txt = "Neoplasia Prev OU", type = "fan")
#nodelabels(text = percentage_labels,node=1:tree$Nnode+Ntip(tree), bg = "black",cex = 0.8, frame = "circle")


# Calculate the range of your data focusing on the majority below 2
small_value_range <- c(min(as.numeric(branches * 100)), 2)
large_value_range <- c(2, max(as.numeric(branches * 100)))

# More breaks in the small value range
small_breaks <- seq(from = small_value_range[1], to = small_value_range[2], length.out = 8)

# Fewer breaks above 2, to ensure these get differentiated colors
large_breaks <- seq(from = large_value_range[1], to = large_value_range[2] + 0.1, length.out = 3)  # Adding 0.1 to ensure coverage

# Combine breaks
breaks <- c(small_breaks, large_breaks[-1])  # Exclude the first of the large breaks to avoid duplication

# Define the colors with emphasis on variety below 2
colors_below_2 <- colorRampPalette(c("blue", "green", "yellow"))(length(small_breaks) - 1)
colors_above_2 <- c("orange", "red")  # Direct assignment for above 2

# Combine colors
colors <- c(colors_below_2, colors_above_2)

# Assign colors based on intervals
color_labels <- colors[findInterval(as.numeric(branches * 100), vec = breaks)]

# Visualization of the tree
obj <- contMap(tree, x, method = "user", anc.states = fitOU, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj, fsize = .8, leg.txt = "Neoplasia Prev OU", type = "fan")
# Applying labels and colors to each node
for (i in 1:length(percentage_labels)) {
  nodelabels(text = sprintf("%.2f%%", branches * 100)[i], node = nodeslist[i], cex = 0.3, frame = "circle", bg = color_labels[i])
}
add.color.bar(leg=80, cols=colors, title="Rate Change from Parent (%)", lims=c(min(breaks), signif(max(branches * 100)), digits = 2), digits=2, prompt=FALSE, lwd=4, outline=TRUE, x=260, y=-230)


par("black")
obj <- contMap(tree, x, method = "user", anc.states = fitOU, plot = TRUE)
obj <- setMap(obj, invert = TRUE,topo.colors(n=20))
plot(obj, fsize = .8, leg.txt = "log(Evolutionary Rate) EB")
labelnodes(text = percentage_labels,node=1:tree$Nnode+Ntip(tree),circle.exp = 0.9, interactive = FALSE)




fitOUdf<-as.data.frame(fitOU)
fitOUdf$nodes<-rownames(fitOUdf)

branchdf$nodes<-rownames(branchdf)
evoPrevNode<-left_join(branchdf, fitOUdf, by = "nodes")
colnames(evoPrevNode)<-c("evo", "node", "prev")
rownames(evoPrevNode)<-evoPrevNode$node
evoPrevNode$prev<-log(evoPrevNode$prev)

prevTip<- as.data.frame(x)
prevTip$tip<-rownames(prevTip)

Evotips<- as.data.frame(tips)
Evotips$tip<-rownames(Evotips)

evoPrevTip<-left_join(prevTip, Evotips, by = "tip")
colnames(evoPrevTip)<-c("prev", "tip", "evo")
rownames(evoPrevTip)<-evoPrevTip$tip
evoPrevTip$prev <- ifelse(evoPrevTip$prev == 0, 0.0001, evoPrevTip$prev)
evoPrevTip$prev<-log(evoPrevTip$prev)


evoPrevNode<-evoPrevNode[,c(3,1),drop=FALSE] 
evoPrevTip<-evoPrevTip[,c(1,3),drop=FALSE] 

# Standardization (Z-score normalization)
standardize <- function(x) {
  return ((x - mean(x)) / sd(x))
}

evoPrevNode$evo <- standardize(evoPrevNode$evo)
evoPrevNode$prev <- standardize(evoPrevNode$prev)

evoPrevTip$evo <- standardize(evoPrevTip$evo)
evoPrevTip$prev <- standardize(evoPrevTip$prev)

# cols<-setNames(rainbow(n =),
#                levels(as.factor(orderData$Orders)))


phylomorphospace(tree, evoPrevTip, A=evoPrevNode, label = "off", node.by.map=TRUE,bty="n", cols = cols,las=1 )

data(sunfish.data)
cols<-setNames(palette()[c(4,2)],
               levels(sunfish.data$feeding.mode))

sunfish.data[,3:2]
