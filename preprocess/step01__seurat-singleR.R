##############################################
#### Install and set the working directory
##############################################

wd <- system("pwd",intern=T)
wd <- paste0(wd,'/temp_result/')
setwd(wd)

#install.packages("remotes")
#remotes::install_github("LTLA/celldex")
library('celldex')


##############################################
#### Prepare the combined reference
##############################################

hpca <- HumanPrimaryCellAtlasData()
bpe <- BlueprintEncodeData()

hpca2 <- hpca
hpca2$label.main <- paste0("HPCA.", hpca2$label.main)
bpe2 <- bpe
bpe2$label.main <- paste0("BPE.", bpe2$label.main)

shared <- intersect(rownames(hpca2), rownames(bpe2))
combined <- cbind(hpca2[shared,], bpe2[shared,])


##############################################
#### Prepare the input file and filtering steps
##############################################

library('Seurat')
library('SingleCellExperiment')
library('SingleR')

f_name <- '../scrna.raw.cnt.tsv'
raw.data <- read.csv(f_name, sep="\t", row.names=1)

canc <- CreateSeuratObject(counts = raw.data, project = "scRNA_cancer", min.cells = 3, min.features = 200)
canc[["percent.mt"]] <- PercentageFeatureSet(canc, pattern = "^MT-")

#### Visualize QC metrics as a violin plot
#VlnPlot(canc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(canc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(canc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

#### Change these cutoff-options depending on your data
canc <- subset(canc, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20 & nCount_RNA < 100000)


##############################################
#### Normalization & Dimensional reduction
##############################################

canc <- NormalizeData(canc, normalization.method = "LogNormalize", scale.factor = 10000)
canc <- FindVariableFeatures(canc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(canc)
canc <- ScaleData(canc, features = all.genes)
canc <- RunPCA(canc, features = VariableFeatures(object = canc))
canc <- FindNeighbors(canc, dims = 1:20)
canc <- FindClusters(canc)#, resolution = 0.5)
canc <- RunTSNE(canc, dims = 1:20)
canc <- RunUMAP(canc, dims = 1:20)

#### Making the output
umap = canc@reductions$umap@cell.embeddings
tsne = canc@reductions$tsne@cell.embeddings
cluster = canc@active.ident
write.table(umap, file="scrna.umap.txt", sep="\t", quote=F)
write.table(tsne, file="scrna.tsne.txt", sep="\t", quote=F)
write.table(cluster, file="scrna.cluster.txt", sep="\t", quote=F)
write.table(all.genes, file="scrna.gene.txt", sep="\t", quote=F)


##############################################
#### SingleR
##############################################

canc.sce <- as.SingleCellExperiment(canc)
result <- SingleR(canc.sce, ref=combined, labels=combined$label.main, assay.type.test=1)
write.table(result, file="scrna.celltype.txt", sep="\t", quote=F)

#### Save the output
canc$pruned.call <- result$pruned.labels
canc$CellType <- result$labels
saveRDS(canc, file = "scrna.seurat-object.rds")


##############################################
#### Making the output for too large data
##############################################

library(Matrix)

writeMM(GetAssayData(object = canc, slot = "data"),file='scrna.seurat-norm.log1p.tsv', sep="\t", quote=F, row.names = T)
writeMM(GetAssayData(object = canc, slot = "counts"),file='scrna.seurat-raw.counts.tsv', sep="\t", quote=F, row.names = T)

sparse_matrix = readMM('./scrna.seurat-norm.log1p.tsv')
extract_file = as.data.frame(sparse_matrix)
rownames(extract_file) = rownames(GetAssayData(object = canc, slot = "data"))
colnames(extract_file) = colnames(GetAssayData(object = canc, slot = "data"))
write.table(extract_file, file='scrna.seurat-norm.log1p.tsv', sep="\t", quote=F, row.names = T)

sparse_matrix_1 = readMM('./scrna.seurat-raw.counts.tsv')
extract_file_1 = as.data.frame(sparse_matrix_1)
rownames(extract_file_1) = rownames(GetAssayData(object = canc, slot = "counts"))
colnames(extract_file_1) = colnames(GetAssayData(object = canc, slot = "counts"))
write.table(extract_file_1, file='scrna.seurat-raw.counts.tsv', sep="\t", quote=F, row.names = T)


####

