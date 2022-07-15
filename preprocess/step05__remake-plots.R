####

wd <- system("pwd",intern=T)
wd <- paste0(wd,'/temp_result/')
setwd(wd)

####

library(Seurat)
library(RColorBrewer)
library(ggplot2)

##

all <- read.csv("scrna.all-annot.t-lab.bi.txt",sep="\t",row.names=1)
png(filename="scrna.plot.tsne.celltype.t-lab.bi.png",
    width = 960, height = 960)
colorCount = length(unique(all$CellType))
getPalette = colorRampPalette(brewer.pal(12,"Paired"))
ggplot(all, aes(tSNE.1, tSNE.2, color=factor(CellType))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values=getPalette(colorCount)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16)
  )
dev.off()

##

all <- read.csv("scrna.all-annot.txt",sep="\t",row.names=1)
png(filename="scrna.plot.tsne.celltype.cluster.png",
    width = 960, height = 960)
colorCount = length(unique(all$Cluster))
getPalette = colorRampPalette(brewer.pal(12,"Paired"))
ggplot(all, aes(tSNE.1, tSNE.2, color=factor(Cluster))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values=getPalette(colorCount)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16)
  )
dev.off()

png(filename="scrna.plot.tsne.celltype.png",
    width = 960, height = 960)
colorCount = length(unique(all$CellType))
getPalette = colorRampPalette(brewer.pal(12,"Paired"))
ggplot(all, aes(tSNE.1, tSNE.2, color=factor(CellType))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values=getPalette(colorCount)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16)
  )
dev.off()


##

all <- read.csv("scrna.all-annot.t-lab.txt",sep="\t",row.names=1)
png(filename="scrna.plot.tsne.celltype.t-lab.png",
    width = 960, height = 960)
colorCount = length(unique(all$CellType))
getPalette = colorRampPalette(brewer.pal(12,"Paired"))
ggplot(all, aes(tSNE.1, tSNE.2, color=factor(CellType))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values=getPalette(colorCount)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16)
  )
dev.off()

##

all <- read.csv("scrna.all-annot.t-lab.txt",sep="\t",row.names=1)
png(filename="scrna.plot.tsne.celltype.batch.png",
    width = 960, height = 960)
colorCount = length(unique(all$Batch))
getPalette = colorRampPalette(brewer.pal(12,"Paired"))
ggplot(all, aes(tSNE.1, tSNE.2, color=factor(Batch))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(size=0.5) +
  theme_classic() +
  scale_color_manual(values=getPalette(colorCount)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16)
  )
dev.off()

####

