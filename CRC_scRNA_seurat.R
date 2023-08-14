setwd() # set work directory
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)


# Load count matrix
cells <- read.csv("GSE146409/GSE146409_UMI_counts_of_filtered_cells.csv", row.names = 1)
# Load metadata
meta <- read.csv("GSE146409/GSE146409_metadata_of_filtered_cells.csv", row.names = 1)

# Create Seurat object
liver <- CreateSeuratObject(counts = cells, meta.data = meta, project = "liver", min.cells = 3, min.features = 200)
# Quality control
liver <- subset(liver, nFeature_RNA > 200 & nFeature_RNA < 10000 & percentMt < 35)
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

# Cluster the cells
liver <- FindNeighbors(liver, dims = 1:15)
liver <- FindClusters(liver, resolution = 1)

# Run non-linear dimensional reduction (tSNE)
liver <- RunTSNE(liver, dims = 1:15)
DimPlot(liver, reduction = "tsne", label = T, label.size = 3)

# Find markers 
markers <- FindAllMarkers(liver, only.pos = TRUE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25,
                          test.use = "wilcox")

top10 <- markers %>% group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) 

## Determine cell types by markers according to The Human Protein Atlas and literature
## Prepare an assignment table 
# Assign cell type to cluster
dd <- read.csv("CRC_cell_type_assign_table.csv")
dd$cluster <- as.factor(dd$cluster)
new.id <- dd$cellType 
names(new.id) <- levels(liver)
liver <- RenameIdents(liver, new.id)
# Add cell type to the metadata as a new column
clst <- data.frame(cluster=liver$seurat_clusters)
clst <- clst %>% left_join(dd, by="cluster")
liver$newIdent <- clst$cellType

# characteristic genes of interest
genes1 <- c("TIMD4","CD14","CSF1R","SPI1")
# Feature plot to show expression pattern
FeaturePlot(liver, features = genes1, reduction = "tsne",
            ncol = 3, cols = c("light grey", "#C11B17"), 
            order = TRUE, keep.scale = "all")

# plot tSNE clusters using ggplot2
md <- liver[[]]
coords <- Embeddings(liver[["tsne"]])
data <- cbind(md, coords)
cols <- c(brewer.pal(7, "Paired"), brewer.pal(8, "Dark2"))
# adjust the order of colors
cols <- c("#E31A1C", "#1F78B4", "#A6761D",  "#B2DF8A", "#FDBF6F",
          "#E7298A", "#E6AB02", "#1B9E77",  "#D95F02", "#A6CEE3",
          "#FB9A99", "#33A02C", "#7570B3",  "#66A61E", "#666666")
ggplot(data, aes(x = tSNE_1, y = tSNE_2, fill = newIdent)) +
  geom_point(shape = 21, size = 1.0, stroke = 0.2, alpha = 0.7) + 
  theme_classic() + theme(legend.position = "None") + 
  xlim(-50,50) + ylim(-50,50) +
  scale_fill_manual(values = cols) 
# plot non-tumor individual
data %>% filter(human == "p6") %>% 
  ggplot(aes(x = tSNE_1, y = tSNE_2, fill = newIdent)) +
  geom_point(shape = 21, size = 1.0, stroke = 0.2, alpha = 0.7) + 
  theme_classic() + theme(legend.position = "None") + 
  xlim(-50,50) + ylim(-50,50) +
  scale_fill_manual(values = cols) 
# plot tumor individuals
data %>% filter(human != "p6") %>% 
  ggplot(aes(x = tSNE_1, y = tSNE_2, fill = newIdent)) +
  geom_point(shape = 21, size = 1.0, stroke = 0.2, alpha = 0.7) + 
  theme_classic() + theme(legend.position = "None") + 
  xlim(-50,50) + ylim(-50,50) +
  scale_fill_manual(values = cols) 


# Violin plot to compare gene expression between KCs and TAMs
genes2 <- c("TIMD4","SIRPA","ID3","CCL3","CCL4","IL18")
# extract counts data
counts <- GetAssayData(object = liver, slot = "data")
gene.use <- match(genes2, rownames(counts))
KC <- WhichCells(object = liver, ident = "KCs")
TAM <- WhichCells(object = liver, ident = "TAMs")
expr.KC <- counts[gene.use, KC]
expr.TAM <- counts[gene.use, TAM]
expr.KC <- as(Class = 'matrix', object = expr.KC)
expr.TAM <- as(Class = 'matrix', object = expr.TAM)
expr.KC <- as.data.frame(t(expr.KC))
expr.TAM <- as.data.frame(t(expr.TAM))
# prepare a data frame for plotting
x <- gather(expr.KC, key = "gene", value = "expression")
x$cell <- rep(rownames(expr.KC), ncol(expr.KC))
x$identity <- rep("KCs", nrow(x))
y <- gather(expr.TAM, key = "gene", value = "expression")
y$cell <- rep(rownames(expr.TAM), ncol(expr.TAM))
y$identity <- rep("TAMs", nrow(y))
expr <- rbind(x, y)
# violin plot using ggplot2 
ggplot(expr, aes(x=identity, y=expression)) + 
  geom_violin(aes(fill=identity)) + 
  stat_summary(fun=mean, geom="point") +
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
        axis.ticks = element_line()) 


# Statistics 
# get the number of cells of each cell type in each individual
cell_num <- data %>% group_by(human, newIdent) %>% summarise(n=n())
cell_num <- spread(cell_num, key = "human", value = "n")

# get the mean expression of genes (normalized)
colMeans(expr.KC)
colMeans(expr.TAM)
# get adjusted p-values between KCs and TAMs 
pv <- FindMarkers(liver, ident.1 = "KCs", ident.2 = "TAMs", features = genes2, 
                  only.pos = FALSE, logfc.threshold = 0, min.pct = 0)


