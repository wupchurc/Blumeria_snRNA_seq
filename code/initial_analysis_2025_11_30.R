# load libraries
library(Seurat)
library(tidyverse)

# Load dataset from 10X CellRanger 
data_dir <- "data/Analysis/cellranger_Control_SP_11_GEX_FL-Z0041/raw_feature_bc_matrix/"
mtx_obj <- Read10X(data.dir = data_dir)
# Initialize seurat object with the raw data
seurat_obj <- CreateSeuratObject(counts = mtx_obj, min.cells = 3, min.features = 200)

# QC -----

# MT
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
        , ncol = 3)
print(p1)
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')
print(p2)

# filter
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data -----
seurat_obj <- NormalizeData(seurat_obj)

# Identify highly variable features -----
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
p3 <- VariableFeaturePlot(seurat_obj)
p4 <- LabelPoints(plot = p3,
                  points = top10,
                  repel = FALSE)
print(p3)
print(p4)

# Scale
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

#Perform linear dimensional reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
# Examine and visualize PCA results a few different ways
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
p5 <- VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
p6 <- DimPlot(seurat_obj, reduction = "pca") + NoLegend()
print(p5)
print(p6)
p7 <- DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
print(p7)
print(ElbowPlot(seurat_obj))

# Cluster the cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
print(head(Idents(seurat_obj), 5))

# Run non-linear dimensional reduction (UMAP/tSNE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p8 <- DimPlot(seurat_obj, reduction = "umap")
print(p8)

# Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurat_obj, ident.1 = 2)
print(head(cluster2.markers, n = 5))
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
seurat_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

