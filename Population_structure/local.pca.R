
library(lostruct)
library(tidyr)
library(ggplot2)
require(data.table)
library(dplyr)
library(stringr)

# Add complete PATH to args 1
# Args 2 corresponds to numerical name of chr (e.g. 1-29, 1A, 4A, Z)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide the file name as an argument.")
}

file_name <- args[1]
CHROM <- args[2]

# Check if the provided file exists
if (!file.exists(file_name)) {
  stop("The provided file does not exist.")
}
setwd("/mnt/DATA/B_JUN_reseq/O_structural_variation/B_SpGen")

# overwrite color discrete to be colorblind friendly

Color <- c("cyan", "red", "gold", "purple", "royalblue1", "purple4", "darkred", "aquamarine2", "orange", "orchid1", "green3", "springgreen4", "yellow", "yellow3", "palegreen", "green", "olivedrab", "gray75")

####data
samples_info <- read.csv("/mnt/DATA/B_JUN_reseq/A_chr_split/B_VCFs/E_JUN_227/chromosome/sample_info.csv"); head(samples_info)
samples <- lostruct::vcf_samples(file_name)

snps<- lostruct::vcf_windower(file_name, size =10000, type = "bp", 
                              samples = lostruct::vcf_samples(file_name))

# extract regions  - Ive made XX windows (number of rows) of 1kb size

regions <- region(snps)(); nrow(regions) ; head (regions)

# Obtain all local PCAs 
# with windowed method, add mc.cores to parallelise
# stores PC1 and PC3 of all local pcas, each containing all individuals
pcs <- eigen_windows(snps, k =2, mc.cores = 18) # can use the option mc.cores = X	

######## MDS of local PCAs

# creates distance matrix of all local PCAs (VERY SLOW)
pcadist <- pc_dist(pcs, npc = 2, mc.cores = 18)
#generates a pairwise distance matrix for mds analysis

#if there are some local PCAs with NAs
na.inds <- is.na(pcadist[,1 ])
na.inds.NAME <- paste0("/mnt/DATA/B_JUN_reseq/O_structural_variation/B_SpGen/na.inds_chr", CHROM)
saveRDS(na.inds, na.inds.NAME)

# next, perform a 2D MDS - slow
# some may have NAs
fit2d <- cmdscale(pcadist[!na.inds,!na.inds], eig=TRUE, k=2)
fit2d.NAME <- paste0("/mnt/DATA/B_JUN_reseq/O_structural_variation/B_SpGen/fit2d_chr", CHROM)
saveRDS(fit2d,  fit2d.NAME)

#Handling MDS: corners
# find extremes
mds_coords <- fit2d$points
mds_corners <- corners(mds_coords, prop = 0.05)

# will always find corners but if very small, this will be a single row, 
# so need to convert to a matrix
if(is.null(dim(mds_corners))){mds_corners <- matrix(mds_corners) %>% t()} 
           
#Make MDS coords into a data.frame
mds <- as_tibble(mds_coords) %>% rename(mds1 = V1, mds2 = V2); head(mds)

# add points vector
mds <- mds %>% mutate(outlier = rep("point", nrow(mds)))
# now add corners
mds$outlier[mds_corners[,1]] <- "corner1"
mds$outlier[mds_corners[,2]] <- "corner2"
mds$outlier[mds_corners[,3]] <- "corner3"

# Add regions and make middle 
mds <- bind_cols(regions[!na.inds,], mds)
mds <- mds %>% mutate(mid = start+(end-start)/2)

# plotting MDS
outlier_colors <- c("corner1" = "#56B4E9", "point" = "gray", "corner3" = "#009E73", "corner2" = "#CC79A7")

# Create a scatter plot with mds1 and mds2, colored by outlier
MDSplot <- ggplot(mds, aes(x = mds1, y = mds2, color = outlier)) +
  geom_point(size = 1) +  # Adjust the size here
  geom_point() +
  scale_color_manual(values = outlier_colors) +
  labs(x = "MDS1", y = "MDS2", color = "Outlier") +
  theme_classic()
MDSplot_name <- paste0("MDS1vsMDS2_chr", CHROM, ".tiff")
# Save the plot as TIFF
ggsave(filename = MDSplot_name, plot = MDSplot, width = 6, height = 5, units = "in", dpi = 300)

MDS1_POS_plot <-  ggplot(mds, aes(mid/10^6, mds1, colour = outlier)) +
  geom_point(size = 1) +
  geom_point() +
  scale_color_manual(values = outlier_colors) +
  labs(x = "midPos (x10⁶)", y = "MDS1", color = "Outlier") +
  theme_classic()
MDS1_POS_plot_name <- paste0("MDS1vsPOS_chr", CHROM, ".tiff")
# Save the plot as TIFF
ggsave(filename = MDS1_POS_plot_name, plot = MDS1_POS_plot, width = 10, height = 5, units = "in", dpi = 300)

MDS2_POS_plot <-  ggplot(mds, aes(mid/10^6, mds2, colour = outlier)) +
  geom_point(size = 1) +
  geom_point() +
  scale_color_manual(values = outlier_colors) +
  labs(x = "midPos (x10⁶)", y = "MDS2", color = "Outlier") +
  theme_classic()
MDS2_POS_plot_name <- paste0("MDS2vsPOS_chr", CHROM, ".tiff")
# Save the plot as TIFF
ggsave(filename = MDS2_POS_plot_name, plot = MDS2_POS_plot, width = 10, height = 5, units = "in", dpi = 300)

##### Corners of MDS
#### PCA per corner
# find the windows - just do separately, it's easier
# Getting SNPs and eigenvectors for corner 1
corner1_snps <- snps(mds_corners[, 1])
# perform eigen window analysis 
corner1_eig <- eigen_windows(corner1_snps, k=2, win =nrow(corner1_snps), mc.cores = 10)

# Getting SNPs and eigenvectors for corner 2
corner2_snps <- snps(mds_corners[, 2])
# perform eigen window analysis 
corner2_eig <- eigen_windows(corner2_snps, k=2, win =nrow(corner2_snps), mc.cores = 10)

# Getting SNPs and eigenvectors for corner 1
corner3_snps <- snps(mds_corners[, 3])
# perform eigen window analysis 
corner3_eig <- eigen_windows(corner3_snps, k=2, win =nrow(corner3_snps), mc.cores = 10)

create_pca <- function(x, samples, sample_info){
  my_pc<- matrix(x[1, 4:ncol(x)], nrow=nrow(samples_info) )
  # add samples and make tibble
  my_pc<- as_tibble(data.frame(ind.vcf = samples, PC1 = my_pc[, 1], PC2 = my_pc[, 2]))
  #combine with sample info to aid plotting
  my_pca <- left_join(my_pc, sample_info, by = "ind.vcf")
  return(my_pca)
}

# Make PCA
corner1_pca <-create_pca(corner1_eig, samples, samples_info)
corner2_pca <-create_pca(corner2_eig, samples, samples_info)
corner3_pca <-create_pca(corner3_eig, samples, samples_info)

PCA_corner1_plot <- ggplot(corner1_pca, aes(PC1, PC2, colour = spp))+
  geom_point(size=3)+
  scale_size_manual(values = c(3, 1)) +
  theme_classic() +
  scale_color_manual(values=Color)+
  ggtitle("Corner 1")
PCA_corner1_plot_name <- paste0("PCA_MDScorner1_chr", CHROM, ".tiff")
# Save the plot as TIFF
ggsave(filename = PCA_corner1_plot_name, plot = PCA_corner1_plot, width = 6, height = 5, units = "in", dpi = 300)

PCA_corner2_plot <- ggplot(corner2_pca, aes(PC1, PC2, colour = spp))+
  geom_point(size=3)+
  scale_size_manual(values = c(3, 1)) +
  theme_classic() +
  scale_color_manual(values=Color)+
  ggtitle("Corner 2")
PCA_corner2_plot_name <- paste0("PCA_MDScorner2_chr", CHROM, ".tiff")
# Save the plot as TIFF
ggsave(filename = PCA_corner2_plot_name, plot = PCA_corner2_plot, width = 6, height = 5, units = "in", dpi = 300)

PCA_corner3_plot <- ggplot(corner3_pca, aes(PC1, PC2, colour = spp))+
  geom_point(size=3)+
  scale_size_manual(values = c(3, 1)) +
  theme_classic() +
  scale_color_manual(values=Color)+
  ggtitle("Corner 3")
PCA_corner3_plot_name <- paste0("PCA_MDScorner3_chr", CHROM, ".tiff")
# Save the plot as TIFF
ggsave(filename = PCA_corner3_plot_name, plot = PCA_corner3_plot, width = 6, height = 5, units = "in", dpi = 300)
