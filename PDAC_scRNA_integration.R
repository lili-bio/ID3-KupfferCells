setwd() # set work directory
library(tidyverse)
library(Seurat)
library(harmony)
library(ggplot2)

# create Seurat object list of cancer patients
samples <- c("P25", "P21", "P27")
seurat_list <- lapply(samples, function(sample){
  cur_data <- Read10X(data.dir = paste0("GSE205013/", sample))
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project=sample
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})

# create Seurat object of healthy individuals
cells <- read.csv("GSE115469/GSE115469_Data.csv", row.names = 1)
health <- CreateSeuratObject(counts = cells, project = "health", min.cells = 3, min.features = 200)
health$SampleID <- health@active.ident

# combine Seurat objects
seurat_list <- append(seurat_list, health)

# merge Seurat objects
merged_seurat <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
# clean up
rm(seurat_list)

# quality control
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat <- subset(merged_seurat, nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# Perform log-normalization and feature selection 
merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
# perform SCT normalization on global object
merged_seurat <- merged_seurat %>% SCTransform()

# Calculate PCs using variable features determined by SCTransform (3000 by default)
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)

# Run Harmony to integrate the patients
harmonized_seurat <- RunHarmony(merged_seurat, 
                                group.by.vars = c("SampleID"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
# run tsne with top PCs
harmonized_seurat <- RunTSNE(harmonized_seurat, reduction='harmony', dims = 1:40)
harmonized_seurat <- FindNeighbors(harmonized_seurat, reduction='harmony')
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 1.8)
DimPlot(harmonized_seurat, reduction = "tsne", group.by = "SampleID")

## The codes for cell type assignment can be found in the CRC script
## KCs were defined as CSF1R+ SPI1+ TIMD4+. TAMs were defined as CSF1R+ SPI1+ TIMD4- TREM2+

## Here we use the pre-saved rds file for convenience. This file is available upon request
harmonized_seurat <- readRDS("PDAC_integration_seurat.rds")

DimPlot(harmonized_seurat, reduction = "tsne", group.by = "newIdent", label = TRUE)

# Feature plot
genes1 <- c("TIMD4","CD14","CSF1R","SPI1")
FeaturePlot(harmonized_seurat, features = genes1, reduction = "tsne",
            ncol = 3, cols = c("light grey", "#C11B17"), 
            order = TRUE, keep.scale = "all")


## violin plot using ggplot
# Extract counts data and prepare a data frame
counts <- GetAssayData(object = harmonized_seurat, slot = "data", assay = "RNA")
genes2 <- c("TIMD4","SIRPA","ID3","CCL3","CCL4","IL18")
gene.use <- match(genes2, rownames(counts))
KC <- WhichCells(object = harmonized_seurat, ident = "KCs")
TAMs <- WhichCells(object = harmonized_seurat, ident = "TAMs")
expr.KC <- counts[gene.use, KC]
expr.TAMs <- counts[gene.use, TAMs]
expr.KC <- as(Class = 'matrix', object = expr.KC)
expr.TAMs <- as(Class = 'matrix', object = expr.TAMs)
expr.KC <- as.data.frame(t(expr.KC))
expr.TAMs<- as.data.frame(t(expr.TAMs))
x <- gather(expr.KC, key = "gene", value = "expression")
x$cell <- rep(rownames(expr.KC), ncol(expr.KC))
x$identity <- rep("KCs", nrow(x))
y <- gather(expr.TAMs, key = "gene", value = "expression")
y$cell <- rep(rownames(expr.TAMs), ncol(expr.TAMs))
y$identity <- rep("TAMs", nrow(y))
expr.comb <- rbind(x,y)

ggplot(expr.comb, aes(x=identity, y=expression)) + 
  geom_violin(aes(fill=identity)) + 
  stat_summary(fun=mean, geom="point", color="black") +
  geom_jitter(size=0.1) +
  scale_fill_manual(values = c("#E7298A", "#7570B3")) +
  labs(y="Expression level") +
  facet_wrap(~gene, scales = "free", ncol = 4) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")) 


# Statistics 
# extract data from seurat object
md <- harmonized_seurat[[]]
# get the number of cells of each cell type in each individual
cell_num <- md %>% group_by(orig.ident, newIdent) %>% summarise(n=n())
cell_num <- spread(cell_num, key = "orig.ident", value = "n")

# get the mean expression of genes (normalized)
colMeans(expr.KC)
colMeans(expr.TAMs)
# get adjusted p-values between KCs and TAMs 
pv <- FindMarkers(harmonized_seurat, ident.1 = "KCs", ident.2 = "TAMs", 
                  features = genes2, assay = "RNA",
                  only.pos = FALSE, logfc.threshold = 0, min.pct = 0)





