
library(ggplot2)
library(data.table)
library(ape)
library(phytools)
library(ggtree)
library(phangorn)  # To read distance matrices



acacti <- fread('/mnt/DATA/B_JUN_reseq/I_phylogeny/G_local_phy/A_ASIP/local_classification.bed') %>%
  setnames(c("chr","start","end","length","cactus","likelihood"))

p <- ggplot(acacti) +
  geom_rect(aes(xmin=start, xmax=end, ymin=-1, ymax=1, fill=cactus)) +
  theme_classic() +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  theme_classic()

print(p)

# Read the distance matrix files and build NJ trees
trees <- lapply(list.files('/mnt/DATA/B_JUN_reseq/I_phylogeny/G_local_phy/A_ASIP/phylip', full.names=TRUE), function(file) {
  dist_matrix <- readDist(file)  # Use phangorn's readDist to read the distance matrix
  tree <- nj(dist_matrix)  # Build a tree using the Neighbor Joining method
  root(tree, outgroup="VUL1", resolve.root=TRUE)  # Root the tree with the specified outgroup
})

# Plot example cactus

ggtree(trees[[1]]) + geom_tiplab() + theme_tree()
ggtree(trees[[2]]) + geom_tiplab() + theme_tree()
ggtree(trees[[3]]) + geom_tiplab() + theme_tree()
ggtree(trees[[4]]) + geom_tiplab() + theme_tree()
ggtree(trees[[5]]) + geom_tiplab() + theme_tree()
ggtree(trees[[6]]) + geom_tiplab() + theme_tree()
ggtree(trees[[7]]) + geom_tiplab() + theme_tree()
ggtree(trees[[8]]) + geom_tiplab() + theme_tree()
ggtree(trees[[9]]) + geom_tiplab() + theme_tree()
ggtree(trees[[10]]) + geom_tiplab() + theme_tree()
ggtree(trees[[11]]) + geom_tiplab() + theme_tree()
ggtree(trees[[12]]) + geom_tiplab() + theme_tree()
ggtree(trees[[13]]) + geom_tiplab() + theme_tree()
ggtree(trees[[14]]) + geom_tiplab() + theme_tree()
ggtree(trees[[15]]) + geom_tiplab() + theme_tree()
