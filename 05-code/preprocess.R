library(Seurat)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(tidyverse)
library(DoubletFinder)

# Load count matrix data from CellRanger ----
# Control
ctrl.11.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_Control_SP_11_GEX_FL-Z0041/filtered_feature_bc_matrix/")
ctrl.12.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_Control_SP_12_GEX_FL-Z0044/filtered_feature_bc_matrix/")
ctrl.14.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_Control_SP_14_GEX_FL-Z0047/filtered_feature_bc_matrix/")
ctrl.15.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_Control_SP_15_GEX_FL-Z0054/filtered_feature_bc_matrix/")

# MCT-Water
water.1.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Water_1_GEX_FL-Z0043/filtered_feature_bc_matrix/")
water.5.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Water_5_GEX_FL-Z0046/filtered_feature_bc_matrix/")
water.6.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Water_6_GEX_FL-Z0053/filtered_feature_bc_matrix/")
water.9.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Water_9_GEX_FL-Z0056/filtered_feature_bc_matrix/")

# MCT-Blumeria
blum.1.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_1_GEX_FL-Z0042/filtered_feature_bc_matrix/")
blum.3.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_3_GEX_FL-Z0045/filtered_feature_bc_matrix/")
blum.4.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_4_GEX_FL-Z0048/filtered_feature_bc_matrix/")
blum.7.mtx <- Read10X(data.dir = "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_7_GEX_FL-Z0055/filtered_feature_bc_matrix/")

# Create Seurat Objects from Count Matrices ----
ctrl.11.seu <- CreateSeuratObject(counts = ctrl.11.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control11")
ctrl.12.seu <- CreateSeuratObject(counts = ctrl.12.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control12")
ctrl.14.seu <- CreateSeuratObject(counts = ctrl.14.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control14")
ctrl.15.seu <- CreateSeuratObject(counts = ctrl.15.mtx, min.cells = 3, 
                                  min.features = 200, project = "Control15")
water.1.seu <- CreateSeuratObject(counts = water.1.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water1")
water.5.seu <- CreateSeuratObject(counts = water.5.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water5")
water.6.seu <- CreateSeuratObject(counts = water.6.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water6")
water.9.seu <- CreateSeuratObject(counts = water.9.mtx, min.cells = 3, 
                                  min.features = 200, project = "Water9")
blum.1.seu <- CreateSeuratObject(counts = blum.1.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria1")
blum.3.seu <- CreateSeuratObject(counts = blum.3.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria3")
blum.4.seu <- CreateSeuratObject(counts = blum.4.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria4")
blum.7.seu <- CreateSeuratObject(counts = blum.7.mtx, min.cells = 3, 
                                 min.features = 200, project = "Blumeria7")

cell_recov_count = list(length(Cells(ctrl.11.seu)), length(Cells(ctrl.12.seu)),
                        length(Cells(ctrl.14.seu)), length(Cells(ctrl.15.seu)),
                        length(Cells(water.1.seu)), length(Cells(water.5.seu)),
                        length(Cells(water.6.seu)), length(Cells(water.9.seu)),
                        length(Cells(blum.1.seu)), length(Cells(blum.3.seu)),
                        length(Cells(blum.4.seu)), length(Cells(blum.7.seu)))
dblt_exp_prcnt <- as.numeric(unlist(cell_recov_count)) / 1000 * 0.8


# Filtering ----
# Check mitochondrial DNA
ctrl.11.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.11.seu, pattern = "^Mt-")
ctrl.12.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.12.seu, pattern = "^Mt-")
ctrl.14.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.14.seu, pattern = "^Mt-")
ctrl.15.seu[['percent.mt']] <- PercentageFeatureSet(ctrl.15.seu, pattern = "^Mt-")
water.1.seu[['percent.mt']] <- PercentageFeatureSet(water.1.seu, pattern = "^Mt-")
water.5.seu[['percent.mt']] <- PercentageFeatureSet(water.5.seu, pattern = "^Mt-")
water.6.seu[['percent.mt']] <- PercentageFeatureSet(water.6.seu, pattern = "^Mt-")
water.9.seu[['percent.mt']] <- PercentageFeatureSet(water.9.seu, pattern = "^Mt-")
blum.1.seu[['percent.mt']] <- PercentageFeatureSet(blum.1.seu, pattern = "^Mt-")
blum.3.seu[['percent.mt']] <- PercentageFeatureSet(blum.3.seu, pattern = "^Mt-")
blum.4.seu[['percent.mt']] <- PercentageFeatureSet(blum.4.seu, pattern = "^Mt-")
blum.7.seu[['percent.mt']] <- PercentageFeatureSet(blum.7.seu, pattern = "^Mt-")

#Violin Plots
VlnPlot(ctrl.11.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(ctrl.11.seu$percent.mt, probs = 0.99)
quantile(ctrl.11.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(ctrl.12.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(ctrl.12.seu$percent.mt, probs = 0.99)
quantile(ctrl.12.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(ctrl.14.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(ctrl.14.seu$percent.mt, probs = 0.99)
quantile(ctrl.14.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(ctrl.15.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(ctrl.15.seu$percent.mt, probs = 0.99)
quantile(ctrl.15.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(water.1.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(water.1.seu$percent.mt, probs = 0.99)
quantile(water.1.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(water.5.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(water.5.seu$percent.mt, probs = 0.99)
quantile(water.5.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(water.6.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(water.6.seu$percent.mt, probs = 0.99)
quantile(water.6.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(water.9.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(water.9.seu$percent.mt, probs = 0.99)
quantile(water.9.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(blum.1.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(blum.1.seu$percent.mt, probs = 0.99)
quantile(blum.1.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(blum.3.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(blum.3.seu$percent.mt, probs = 0.99)
quantile(blum.3.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(blum.4.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(blum.4.seu$percent.mt, probs = 0.99)
quantile(blum.4.seu$nFeature_RNA, c(0.05, 0.99))

VlnPlot(blum.7.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(blum.7.seu$percent.mt, probs = 0.99)
quantile(blum.7.seu$nFeature_RNA, c(0.05, 0.99))

# Subsetting/Filtering
ctrl.11.seu <- subset(ctrl.11.seu, subset = nFeature_RNA < 4103.06 & percent.mt < 2.882178)
ctrl.12.seu <- subset(ctrl.12.seu, subset = nFeature_RNA < 4040.52 & percent.mt < 1.348935)
ctrl.14.seu <- subset(ctrl.14.seu, subset = nFeature_RNA < 4709.74 & percent.mt < 1.20078)
ctrl.15.seu <- subset(ctrl.15.seu, subset = nFeature_RNA < 3922.04 & percent.mt < 2.253213)
water.1.seu <- subset(water.1.seu, subset = nFeature_RNA < 4732.02 & percent.mt < 5.821339)
water.5.seu <- subset(water.5.seu, subset = nFeature_RNA < 4448 & percent.mt < 5.646684)
water.6.seu <- subset(water.6.seu, subset = nFeature_RNA < 4719.97 & percent.mt < 3.978664)
water.9.seu <- subset(water.9.seu, subset = nFeature_RNA < 4259.4 & percent.mt < 2.022335 )
blum.1.seu <- subset(blum.1.seu, subset = nFeature_RNA < 3931.39 & percent.mt < 2.343205)
blum.3.seu <- subset(blum.3.seu, subset = nFeature_RNA < 3987.5 & percent.mt < 3.920576)
blum.4.seu <- subset(blum.4.seu, subset = nFeature_RNA < 4846.44 & percent.mt < 5.313926)
blum.7.seu <- subset(blum.7.seu, subset = nFeature_RNA < 4489.29 & percent.mt < 3.517909)

# Normalize ----
ctrl.11.seu <- NormalizeData(ctrl.11.seu)
ctrl.12.seu <- NormalizeData(ctrl.12.seu)
ctrl.14.seu <- NormalizeData(ctrl.14.seu)
ctrl.15.seu <- NormalizeData(ctrl.15.seu)
water.1.seu <- NormalizeData(water.1.seu)
water.5.seu <- NormalizeData(water.5.seu)
water.6.seu <- NormalizeData(water.6.seu)
water.9.seu <- NormalizeData(water.9.seu)
blum.1.seu <- NormalizeData(blum.1.seu)
blum.3.seu <- NormalizeData(blum.3.seu)
blum.4.seu <- NormalizeData(blum.4.seu)
blum.7.seu <- NormalizeData(blum.7.seu)

# Find Variable Features ----
ctrl.11.seu <- FindVariableFeatures(ctrl.11.seu)
ctrl.12.seu <- FindVariableFeatures(ctrl.12.seu)
ctrl.14.seu <- FindVariableFeatures(ctrl.14.seu)
ctrl.15.seu <- FindVariableFeatures(ctrl.15.seu)
water.1.seu <- FindVariableFeatures(water.1.seu)
water.5.seu <- FindVariableFeatures(water.5.seu)
water.6.seu <- FindVariableFeatures(water.6.seu)
water.9.seu <- FindVariableFeatures(water.9.seu)
blum.1.seu <- FindVariableFeatures(blum.1.seu)
blum.3.seu <- FindVariableFeatures(blum.3.seu)
blum.4.seu <- FindVariableFeatures(blum.4.seu)
blum.7.seu <- FindVariableFeatures(blum.7.seu)

# Scale Data ----
ctrl.11.seu <- ScaleData(ctrl.11.seu, feature = rownames(ctrl.11.seu))
ctrl.12.seu <- ScaleData(ctrl.12.seu, feature = rownames(ctrl.12.seu))
ctrl.14.seu <- ScaleData(ctrl.14.seu, feature = rownames(ctrl.14.seu))
ctrl.15.seu <- ScaleData(ctrl.15.seu, feature = rownames(ctrl.15.seu))
water.1.seu <- ScaleData(water.1.seu, feature = rownames(water.1.seu))
water.5.seu <- ScaleData(water.5.seu, feature = rownames(water.5.seu))
water.6.seu <- ScaleData(water.6.seu, feature = rownames(water.6.seu))
water.9.seu <- ScaleData(water.9.seu, feature = rownames(water.9.seu))
blum.1.seu <- ScaleData(blum.1.seu, feature = rownames(blum.1.seu))
blum.3.seu <- ScaleData(blum.3.seu, feature = rownames(blum.3.seu))
blum.4.seu <- ScaleData(blum.4.seu, feature = rownames(blum.4.seu))
blum.7.seu <- ScaleData(blum.7.seu, feature = rownames(blum.7.seu))

# PCA ----
ctrl.11.seu <- RunPCA(ctrl.11.seu, features = VariableFeatures(object = ctrl.11.seu))
ctrl.12.seu <- RunPCA(ctrl.12.seu, features = VariableFeatures(object = ctrl.12.seu))
ctrl.14.seu <- RunPCA(ctrl.14.seu, features = VariableFeatures(object = ctrl.14.seu))
ctrl.15.seu <- RunPCA(ctrl.15.seu, features = VariableFeatures(object = ctrl.15.seu))
water.1.seu <- RunPCA(water.1.seu, features = VariableFeatures(object = water.1.seu))
water.5.seu <- RunPCA(water.5.seu, features = VariableFeatures(object = water.5.seu))
water.6.seu <- RunPCA(water.6.seu, features = VariableFeatures(object = water.6.seu))
water.9.seu <- RunPCA(water.9.seu, features = VariableFeatures(object = water.9.seu))
blum.1.seu <- RunPCA(blum.1.seu, features = VariableFeatures(object = blum.1.seu))
blum.3.seu <- RunPCA(blum.3.seu, features = VariableFeatures(object = blum.3.seu))
blum.4.seu <- RunPCA(blum.4.seu, features = VariableFeatures(object = blum.4.seu))
blum.7.seu <- RunPCA(blum.7.seu, features = VariableFeatures(object = blum.7.seu))

# Elbow Plots ----
ElbowPlot(ctrl.11.seu)
ElbowPlot(ctrl.12.seu)
ElbowPlot(ctrl.14.seu)
ElbowPlot(ctrl.15.seu)
ElbowPlot(water.1.seu)
ElbowPlot(water.5.seu)
ElbowPlot(water.6.seu)
ElbowPlot(water.9.seu)
ElbowPlot(blum.1.seu)
ElbowPlot(blum.3.seu)
ElbowPlot(blum.4.seu)
ElbowPlot(blum.7.seu)

# Cluster ----
ctrl.11.seu <- FindNeighbors(ctrl.11.seu, dims  = 1:20)
ctrl.11.seu <- FindClusters(ctrl.11.seu)
ctrl.11.seu <- RunUMAP(ctrl.11.seu, dims = 1:20)
ctrl.12.seu <- FindNeighbors(ctrl.12.seu, dims  = 1:20)
ctrl.12.seu <- FindClusters(ctrl.12.seu)
ctrl.12.seu <- RunUMAP(ctrl.12.seu, dims = 1:20)
ctrl.14.seu <- FindNeighbors(ctrl.14.seu, dims  = 1:20)
ctrl.14.seu <- FindClusters(ctrl.14.seu)
ctrl.14.seu <- RunUMAP(ctrl.14.seu, dims = 1:20)
ctrl.15.seu <- FindNeighbors(ctrl.15.seu, dims  = 1:20)
ctrl.15.seu <- FindClusters(ctrl.15.seu)
ctrl.15.seu <- RunUMAP(ctrl.15.seu, dims = 1:20)

water.1.seu <- FindNeighbors(water.1.seu, dims  = 1:20)
water.1.seu <- FindClusters(water.1.seu)
water.1.seu <- RunUMAP(water.1.seu, dims = 1:20)
water.5.seu <- FindNeighbors(water.5.seu, dims  = 1:20)
water.5.seu <- FindClusters(water.5.seu)
water.5.seu <- RunUMAP(water.5.seu, dims = 1:20)
water.6.seu <- FindNeighbors(water.6.seu, dims  = 1:20)
water.6.seu <- FindClusters(water.6.seu)
water.6.seu <- RunUMAP(water.6.seu, dims = 1:20)
water.9.seu <- FindNeighbors(water.9.seu, dims  = 1:20)
water.9.seu <- FindClusters(water.9.seu)
water.9.seu <- RunUMAP(water.9.seu, dims = 1:20)

blum.1.seu <- FindNeighbors(blum.1.seu, dims  = 1:20)
blum.1.seu <- FindClusters(blum.1.seu)
blum.1.seu <- RunUMAP(blum.1.seu, dims = 1:20)
blum.3.seu <- FindNeighbors(blum.3.seu, dims  = 1:20)
blum.3.seu <- FindClusters(blum.3.seu)
blum.3.seu <- RunUMAP(blum.3.seu, dims = 1:20)
blum.4.seu <- FindNeighbors(blum.4.seu, dims  = 1:20)
blum.4.seu <- FindClusters(blum.4.seu)
blum.4.seu <- RunUMAP(blum.4.seu, dims = 1:20)
blum.7.seu <- FindNeighbors(blum.7.seu, dims  = 1:20)
blum.7.seu <- FindClusters(blum.7.seu)
blum.7.seu <- RunUMAP(blum.7.seu, dims = 1:20)

# Doublets ----

# Control11
sweep.res <- paramSweep(ctrl.11.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.11.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.15*length(Cells(ctrl.11.seu))) #Assuming 15% doublet formation rate based on .8% doublets per 1000 cells and 16622 cells recovered
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.11.seu <- doubletFinder(ctrl.11.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.11.seu@meta.data), value = TRUE)

DimPlot(ctrl.11.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.11.seu@meta.data[[df_col]])

#subset data to keep only singlets
# ctrl.11.seu <- subset(ctrl.11.seu, subset = ctrl.11.seu@meta.data[[df_col]] == 'Singlet')
ctrl.11.seu <- ctrl.11.seu[, ctrl.11.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.11.seu@meta.data[[df_col]])

# Control12
sweep.res <- paramSweep(ctrl.12.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.12.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.17*length(Cells(ctrl.12.seu))) #Assuming 17% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.12.seu <- doubletFinder(ctrl.12.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.12.seu@meta.data), value = TRUE)

DimPlot(ctrl.12.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.12.seu@meta.data[[df_col]])

#subset data to keep only singlets
ctrl.12.seu <- ctrl.12.seu[, ctrl.12.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.12.seu@meta.data[[df_col]])

# Control14
sweep.res <- paramSweep(ctrl.14.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.14.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.19*length(Cells(ctrl.14.seu))) #Assuming 19% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.14.seu <- doubletFinder(ctrl.14.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.14.seu@meta.data), value = TRUE)

DimPlot(ctrl.14.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.14.seu@meta.data[[df_col]])

#subset data to keep only singlets
ctrl.14.seu <- ctrl.14.seu[, ctrl.14.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.14.seu@meta.data[[df_col]])

# Control15
sweep.res <- paramSweep(ctrl.15.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- ctrl.15.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.20*length(Cells(ctrl.15.seu))) #Assuming 20% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

ctrl.15.seu <- doubletFinder(ctrl.15.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(ctrl.15.seu@meta.data), value = TRUE)

DimPlot(ctrl.15.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(ctrl.15.seu@meta.data[[df_col]])

#subset data to keep only singlets
ctrl.15.seu <- ctrl.15.seu[, ctrl.15.seu@meta.data[[df_col]] == "Singlet"]

table(ctrl.15.seu@meta.data[[df_col]])

# Water1
sweep.res <- paramSweep(water.1.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.1.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.15*length(Cells(water.1.seu))) #Assuming 15% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.1.seu <- doubletFinder(water.1.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.1.seu@meta.data), value = TRUE)

DimPlot(water.1.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.1.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.1.seu <- water.1.seu[, water.1.seu@meta.data[[df_col]] == "Singlet"]

table(water.1.seu@meta.data[[df_col]])

# Water5
sweep.res <- paramSweep(water.5.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.5.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.19*length(Cells(water.5.seu))) #Assuming 19% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.5.seu <- doubletFinder(water.5.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.5.seu@meta.data), value = TRUE)

DimPlot(water.5.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.5.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.5.seu <- water.5.seu[, water.5.seu@meta.data[[df_col]] == "Singlet"]

table(water.5.seu@meta.data[[df_col]])

# Water6
sweep.res <- paramSweep(water.6.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.6.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.18*length(Cells(water.6.seu))) #Assuming 18% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.6.seu <- doubletFinder(water.6.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.6.seu@meta.data), value = TRUE)

DimPlot(water.6.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.6.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.6.seu <- water.6.seu[, water.6.seu@meta.data[[df_col]] == "Singlet"]

table(water.6.seu@meta.data[[df_col]])

# Water9
sweep.res <- paramSweep(water.9.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- water.9.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.16*length(Cells(water.9.seu))) #Assuming 16% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

water.9.seu <- doubletFinder(water.9.seu, PCs = 1:20, pN = 0.25, pK = pK,
                             nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(water.9.seu@meta.data), value = TRUE)

DimPlot(water.9.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(water.9.seu@meta.data[[df_col]])

#subset data to keep only singlets
water.9.seu <- water.9.seu[, water.9.seu@meta.data[[df_col]] == "Singlet"]

table(water.9.seu@meta.data[[df_col]])

# Blumeria1
sweep.res <- paramSweep(blum.1.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.1.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.17*length(Cells(blum.1.seu))) #Assuming 17% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.1.seu <- doubletFinder(blum.1.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.1.seu@meta.data), value = TRUE)

DimPlot(blum.1.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.1.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.1.seu <- blum.1.seu[, blum.1.seu@meta.data[[df_col]] == "Singlet"]

table(blum.1.seu@meta.data[[df_col]])

# Blumeria3
sweep.res <- paramSweep(blum.3.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.3.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.18*length(Cells(blum.3.seu))) #Assuming 18% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.3.seu <- doubletFinder(blum.3.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.3.seu@meta.data), value = TRUE)

DimPlot(blum.3.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.3.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.3.seu <- blum.3.seu[, blum.3.seu@meta.data[[df_col]] == "Singlet"]

table(blum.3.seu@meta.data[[df_col]])

# Blumeria4
sweep.res <- paramSweep(blum.4.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.4.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.16*length(Cells(blum.4.seu))) #Assuming 16% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.4.seu <- doubletFinder(blum.4.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.4.seu@meta.data), value = TRUE)

DimPlot(blum.4.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.4.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.4.seu <- blum.4.seu[, blum.4.seu@meta.data[[df_col]] == "Singlet"]

table(blum.4.seu@meta.data[[df_col]])

# Blumeria7
sweep.res <- paramSweep(blum.7.seu, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group=1)) + geom_point() + geom_line()

pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- blum.7.seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.20*length(Cells(blum.7.seu))) #Assuming 20% doublet formation rate based on .8% doublets per 1000 cells
nExp_poi.adj <- round(nExp_poi*(1- homotypic.prop))

blum.7.seu <- doubletFinder(blum.7.seu, PCs = 1:20, pN = 0.25, pK = pK,
                            nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

df_col <- grep("^DF\\.classifications_", colnames(blum.7.seu@meta.data), value = TRUE)

DimPlot(blum.7.seu, reduction = 'umap', group.by = df_col, raster = FALSE)
table(blum.7.seu@meta.data[[df_col]])

#subset data to keep only singlets
blum.7.seu <- blum.7.seu[, blum.7.seu@meta.data[[df_col]] == "Singlet"]

table(blum.7.seu@meta.data[[df_col]])

# Add metadata ----

ctrl.11.seu$sample <- "Control11"
ctrl.11.seu$condition <- "Control"

ctrl.12.seu$sample <- "Control12"
ctrl.12.seu$condition <- "Control"

ctrl.14.seu$sample <- "Control14"
ctrl.14.seu$condition <- "Control"

ctrl.15.seu$sample <- "Control15"
ctrl.15.seu$condition <- "Control"

water.1.seu$sample <- "Water1"
water.1.seu$condition <- "MCT-Water"

water.5.seu$sample <- "Water5"
water.5.seu$condition <- "MCT-Water"

water.6.seu$sample <- "Water6"
water.6.seu$condition <- "MCT-Water"

water.9.seu$sample <- "Water9"
water.9.seu$condition <- "MCT-Water"

blum.1.seu$sample <- "Blumeria1"
blum.1.seu$condition <- "MCT-Blumeria"

blum.3.seu$sample <- "Blumeria3"
blum.3.seu$condition <- "MCT-Blumeria"

blum.4.seu$sample <- "Blumeria4"
blum.4.seu$condition <- "MCT-Blumeria"

blum.7.seu$sample <- "Blumeria7"
blum.7.seu$condition <- "MCT-Blumeria"

# Combine ----
seu_merge <- merge(ctrl.11.seu, y = list(ctrl.12.seu, ctrl.14.seu, ctrl.15.seu,
                                         water.1.seu, water.5.seu, water.6.seu, water.9.seu,
                                         blum.1.seu, blum.3.seu, blum.4.seu, blum.7.seu),
                   add.cell.ids = c("Control11","Control12", "Control14", "Control15",
                                    "Water1", "Water5", "Water6", "Water9",
                                    "Blumeria1", "Blumeria3", "Blumeria4", "Blumeria7"),
                   project = ("All.Samples"))

saveRDS(seu_merge, file = "03-analysis_scratch/seu_merge_filtered.rds")
# seu_merge <- readRDS("03-analysis_scratch/seu_merge.rds")

seu_merge <- NormalizeData(seu_merge)
seu_merge <- FindVariableFeatures(seu_merge)
seu_merge <- ScaleData(seu_merge)
seu_merge <- RunPCA(seu_merge)

seu_merge <- IntegrateLayers(
  object = seu_merge,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony"
)

# saveRDS(seu_merge, file = "03-analysis_scratch/seu_integrated.rds")
# seu_merge <- readRDS("03-analysis_scratch/seu_integrated.rds")

seu_merge <- FindNeighbors(seu_merge, reduction = "harmony", dims = 1:20)
seu_merge <- FindClusters(seu_merge, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
seu_merge <- RunUMAP(seu_merge, reduction = "harmony", dims = 1:20)

seu_merge$condition <- factor(
  seu_merge$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)
DimPlot(
  seu_merge,
  reduction = "umap",
  group.by = "condition",
  split.by = "condition",
  label = FALSE
)

DimPlot(
  seu_merge,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  split.by = "condition",
  label = TRUE
)
