setwd() # set work directory
library(Seurat)
library(patchwork)
library(Matrix)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)


# Load counts matrix
cells <- read.csv("GSE146409/GSE146409_UMI_counts_of_filtered_cells.csv", row.names = 1)
# Load metadata
meta <- read.csv("GSE146409/GSE146409_metadata_of_filtered_cells.csv", row.names = 1)

# Create Seurat object
liver <- CreateSeuratObject(counts = cells, meta.data = meta, project = "liver", min.cells = 3, min.features = 200)
# Quality control
liver <- subset(liver, nFeature_RNA > 200 & nFeature_RNA < 3000 & percentMt < 35)
# Subset data to remove p1 and p2
liver <- subset(liver, subset = human != "p1" & human != "p2" )

# Normalize the data
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable genes for PCA 
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)
# Scale the data
all.genes <- rownames(cells)
liver <- ScaleData(liver, features = all.genes)

# Perform linear dimensional reduction
liver <- RunPCA(liver, features = VariableFeatures(object = liver))
DimPlot(liver, reduction = "pca")

# Examine how many components to be included
liver <- JackStraw(liver, num.replicate = 100)
liver <- ScoreJackStraw(liver, dims = 1:20)
JackStrawPlot(liver, dims = 1:20)

# Cluster the cells
liver <- FindNeighbors(liver, dims = 1:15)
liver <- FindClusters(liver, resolution = 1)

# Run non-linear dimensional reduction (tSNE)
liver <- RunTSNE(liver, dims = 1:15)
DimPlot(liver, reduction = "tsne", label = T, label.size = 5)

# Find markers 
markers <- FindAllMarkers(liver, only.pos = TRUE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25,
                          test.use = "wilcox")

top10 <- markers %>% group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) 
DoHeatmap(liver, features = top10$gene, size = 3) + NoLegend()

## Determine cell types by markers according to references
## Prepare an assignment table
# Assign cell type to cluster
dd <- read.csv("cluster_rename.csv")
dd$cluster <- as.factor(dd$cluster)
new.id <- dd$cellType 
names(new.id) <- levels(liver)
liver <- RenameIdents(liver, new.id)
# Add cell type to the metadata as well
clst <- data.frame(cluster=liver$seurat_clusters)
clst <- clst %>% left_join(dd, by="cluster")
liver$newIdent <- clst$cellType

# characteristic genes of interest
genes <- c("TIMD4","SIRPA","ID3","CCL3","CCL4","IL18")

# Feature plot to show expression pattern
FeaturePlot(liver, features = genes, reduction = "tsne",
            ncol = 3, cols = c("light grey", "#C11B17"), 
            order = TRUE, keep.scale = "all")

# Violin plot to compare gene expression between macrophage populations
VlnPlot(liver, idents = c("KCs", "TAMs"), features = genes) +
  stat_summary(fun=mean, geom="point")

# Statistics 
# get the mean expression value of genes
ave <- AverageExpression(liver, features = genes)
# get adjusted p-value 
pv <- FindMarkers(liver, ident.1 = "KCs", ident.2 = "TAMs", features = genes, 
                  only.pos = FALSE, logfc.threshold = 0, min.pct = 0)

# Aesthetics
# tSNE plot with ggplot
md <- liver[[]]
coords <- Embeddings(liver[["tsne"]])
data <- cbind(md, coords)
cols <- c(brewer.pal(8, "Dark2"), brewer.pal(7, "Paired"))
ggplot(data, aes(x = tSNE_1, y = tSNE_2, fill = newIdent)) +
  geom_point(shape = 21, size = 1.0, stroke = 0.2, alpha = 0.7) + 
  theme_classic() + theme(legend.position = "None") + 
  xlim(-50,50) + ylim(-50,50) +
  scale_fill_manual(values = cols) 


# Extract counts data for violin plot
counts <- GetAssayData(object = liver, slot = "data")
gene.use <- match(genes, rownames(counts))
KC <- WhichCells(object = liver3, ident = "KCs")
TAM <- WhichCells(object = liver3, ident = "TAMs")
expr.KC <- counts[gene.use, KC]
expr.TAM <- counts[gene.use, TAM]

expr.KC <- as(Class = 'matrix', object = expr.KC)
expr.TAM <- as(Class = 'matrix', object = expr.TAM)
expr.KC <- as.data.frame(t(expr.KC))
expr.TAM <- as.data.frame(t(expr.TAM))

x <- gather(expr.KC, key = "gene", value = "expression")
x$cell <- rep(rownames(expr.KC), ncol(expr.KC))
x$identity <- rep("KCs", nrow(x))
y <- gather(expr.TAM, key = "gene", value = "expression")
y$cell <- rep(rownames(expr.TAM), ncol(expr.TAM))
y$identity <- rep("TAMs", nrow(y))
expr <- rbind(x, y)

# Violin plot using ggplot 
ggplot(expr, aes(x=identity, y=expression)) + 
  geom_violin(aes(fill=identity)) + 
  stat_summary(fun=mean, geom="point") +
  geom_jitter(size=0.1) +
  scale_fill_manual(values = cols[c(4,3)]) +
  labs(y="Expression level") +
  facet_wrap(~gene, scales = "free", ncol = 4) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line()) 



