#/bin/R

setwd("/mnt/DATA/B_JUN_reseq/C_PopAnalysis/E_plink/G_225/JUN_225.auto.neutral.u")

library(ggplot2)
library(tidyverse)

pca <-  read.table("plink.eigenvec")
eigenval <- scan("./plink.eigenval")

pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

spp <- rep(NA, length(pca$ind))
spp[grep("AIK", pca$ind)] <- "AIK"
spp[grep("ALT", pca$ind)] <- "ALT"
spp[grep("BAI", pca$ind)] <- "BAI"
spp[grep("CAN", pca$ind)] <- "CAN"
spp[grep("CAR", pca$ind)] <- "CAR"
spp[grep("DOE", pca$ind)] <- "DOE"
spp[grep("DOW", pca$ind)] <- "DOW"
spp[grep("FUL", pca$ind)] <- "FUL"
spp[grep("HYE", pca$ind)] <- "HYE"
spp[grep("INS", pca$ind)] <- "INS"
spp[grep("MEA", pca$ind)] <- "MEA"
spp[grep("MON", pca$ind)] <- "MON"
spp[grep("ORE", pca$ind)] <- "ORE"
spp[grep("PAL", pca$ind)] <- "PAL"
spp[grep("PHA", pca$ind)] <- "PHA"
spp[grep("PIN", pca$ind)] <- "PIN"
spp[grep("TOW", pca$ind)] <- "TOW"
spp[grep("THU", pca$ind)] <- "THU"
spp[grep("VUL", pca$ind)] <- "VUL"

ALL <- as_tibble(data.frame(pca, spp))

# first convert to percentage variance explained

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()


# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
#b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = area)) + geom_point(size = 2)
#b <- b + scale_colour_manual(values = c("gray", "pink", "red", "gold", "green", "blue", "skyblue"))
#b <- b + scale_shape_manual(values = c( 16, 17))
#b <- b + coord_equal() + theme_light()
#b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

pca %>% ggplot(aes(PC1, PC2, col = spp))+
  geom_point(alpha=1, size=5)+
  scale_colour_manual(values = c("cyan", "red", "gold", "orchid1", "skyblue2", "purple4", "purple", "darkred", "blue", "orange", "pink", "forestgreen", "palegreen1", "yellow", "yellow3", "springgreen2", "darkseagreen", "olivedrab3", "gray75"))+
  #scale_shape_manual(values = c( 17, 8, 15, 16, 0, 4, 2, 3))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),         # No major gridlines
    panel.grid.minor = element_blank(),         # No minor gridlines
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  theme(legend.position = "bottom", aspect.ratio = 5/7, legend.title = element_text(size=12))
ggsave("PCA_PC1-2.tiff", dpi = 300, device = "tiff")


pca %>% ggplot(aes(PC3, PC4, col = spp))+
  geom_point(alpha=1, size=5)+
  scale_colour_manual(values = c("cyan", "red", "gold", "orchid1", "skyblue2", "purple4", "purple", "darkred", "blue", "orange", "pink", "forestgreen", "palegreen1", "yellow", "yellow3", "springgreen2", "darkseagreen", "olivedrab3", "gray75"))+
  #scale_shape_manual(values = c( 17, 8, 15, 16, 0, 4, 2, 3))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),         # No major gridlines
    panel.grid.minor = element_blank(),         # No minor gridlines
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")) +
  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)"))+
  theme(legend.position = "bottom", aspect.ratio = 5/7, legend.title = element_text(size=12))
ggsave("PCA_PC2-3.tiff", dpi = 300, device = "tiff")


pca %>% ggplot(aes(PC5, PC6, col = spp))+
  geom_point(alpha=1, size=5)+
  scale_colour_manual(values = c("cyan", "red", "gold", "orchid1", "skyblue2", "purple4", "purple", "darkred", "blue", "orange", "pink", "forestgreen", "palegreen1", "yellow", "yellow3", "springgreen2", "darkseagreen", "olivedrab3", "gray75"))+
  #scale_shape_manual(values = c( 17, 8, 15, 16, 0, 4, 2, 3))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),         # No major gridlines
    panel.grid.minor = element_blank(),         # No minor gridlines
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")) +
  xlab(paste0("PC5 (", signif(pve$pve[5], 3), "%)")) + ylab(paste0("PC6 (", signif(pve$pve[6], 3), "%)"))+
  theme(legend.position = "bottom", aspect.ratio = 5/7, legend.title = element_text(size=12))
ggsave("PCA_PC5-6.tiff", dpi = 300, device = "tiff")


pca %>% ggplot(aes(PC7, PC8, col = spp))+
  geom_point(alpha=1, size=5)+
  scale_colour_manual(values = c("cyan", "red", "gold", "orchid1", "skyblue2", "purple4", "purple", "darkred", "blue", "orange", "pink", "forestgreen", "palegreen1", "yellow", "yellow3", "springgreen2", "darkseagreen", "olivedrab3", "gray75"))+
  #scale_shape_manual(values = c( 17, 8, 15, 16, 0, 4, 2, 3))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),         # No major gridlines
    panel.grid.minor = element_blank(),         # No minor gridlines
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")) +
  xlab(paste0("PC7 (", signif(pve$pve[7], 3), "%)")) + ylab(paste0("PC8 (", signif(pve$pve[8], 3), "%)"))+
  theme(legend.position = "bottom", aspect.ratio = 5/7, legend.title = element_text(size=12))
ggsave("PCA_PC7-8.tiff", dpi = 300, device = "tiff")
