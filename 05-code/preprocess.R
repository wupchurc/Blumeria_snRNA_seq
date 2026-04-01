# load libraries ----
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(presto)
library(patchwork)
library(scCustomize)
library(DESeq2)
library(SummarizedExperiment)
library(RColorBrewer)
library(circlize)
library(tidyr)
library(dplyr)

# Load control dataset from 10X CellRanger ----
data_dirs <- c("01-raw_data_scratch/Analysis/cellranger_Control_SP_11_GEX_FL-Z0041/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_Control_SP_12_GEX_FL-Z0044/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_Control_SP_14_GEX_FL-Z0047/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_Control_SP_15_GEX_FL-Z0054/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Water_1_GEX_FL-Z0043/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Water_5_GEX_FL-Z0046/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Water_6_GEX_FL-Z0053/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Water_9_GEX_FL-Z0056/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_1_GEX_FL-Z0042/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_3_GEX_FL-Z0045/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_4_GEX_FL-Z0048/filtered_feature_bc_matrix/",
               "01-raw_data_scratch/Analysis/cellranger_MCT_Blumeria_7_GEX_FL-Z0055/filtered_feature_bc_matrix/")
sample_ids <- c("Control11", "Control12", "Control14", "Control15",
                "Water1", "Water5", "Water6", "Water9",
                "Blumeria1", "Blumeria3", "Blumeria4", "Blumeria7")

# Create seurat objects ----
obj_list <- list()
for (i in seq_along(data_dirs)) {
  seurat_obj <- Read10X(data.dir = data_dirs[i])
  obj_list[[i]] <- CreateSeuratObject(counts = seurat_obj, min.cells = 3, 
                                      min.features = 200, project = sample_ids[i])
}

# Add condition label
for (i in seq_along(obj_list)) {
  sid <- sample_ids[i]
  
  if (grepl("^Control",  sid)) obj_list[[i]]$condition <- "Control"
  if (grepl("^Water",    sid)) obj_list[[i]]$condition <- "MCT-Water"
  if (grepl("^Blumeria", sid)) obj_list[[i]]$condition <- "MCT-Blumeria"
}
# Doublet rate estimation----
cell.recov.counts <- lapply(obj_list, ncol)
cells_ref  <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
rate_ref   <- c(0.4,0.8,1.6,2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)   # in percent
# linear interpolation; extrapolate = TRUE if some samples outside range
dbl.rate <- approx(x = cells_ref,
                   y = rate_ref,
                   xout = cell.recov.counts,
                   method = "linear",
                   rule = 2)$y
names(dbl.rate) <- names(cell.recov.counts)
dbl.rate

# Remove Doublets ----
doubletFinder_obj_list <- list()
for (i in seq_along(obj_list)) {
  seu <- obj_list[[i]]
  
  # QC
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^Mt-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  # Preprocessing (required for DoubletFinder)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  # print(ElbowPlot(seu))
  
  # DoubletFinder pK identification (no t-SNE needed)
  sweep.res <- paramSweep(seu, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # select pK value that corresponds to max bcmvn
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  # double estimation
  nExp <- round(cell.recov.counts[[i]] * (dbl.rate[[i]] / 100))
  # run doubletFinder
  seu <- doubletFinder(seu, PCs = 1:20, pN = 0.25, pK = pK, 
                          nExp = nExp, reuse.pANN = NULL, sct = FALSE)

  doubletFinder_obj_list[[i]] <- seu  
}

saveRDS(doubletFinder_obj_list, file = "doubletFinder_obj.rds")
# doubletFinder_obj_list <- readRDS("doubletFinder_obj.rds")

# Filter ----
seu_singlets_list <- vector("list", length(doubletFinder_obj_list))

for (i in seq_along(doubletFinder_obj_list)) {
  seu <- doubletFinder_obj_list[[i]]
  
  df_cols <- grep("^DF.classifications", colnames(seu@meta.data), value = TRUE)
  if (length(df_cols) != 1) {
    stop("Expected exactly one DF.classifications column, found: ",
         paste(df_cols, collapse = ", "))
  }
  df_col <- df_cols[1]
  
  # logical index of singlets
  singlet_idx <- seu@meta.data[[df_col]] == "Singlet"
  seu_singlets <- seu[, singlet_idx]
  
  seu_singlets_list[[i]] <- seu_singlets
}

# assume length(sample_ids) == length(seu_singlets_list)
for (i in seq_along(seu_singlets_list)) {
  seu <- seu_singlets_list[[i]]
  prefix <- sample_ids[i]
  colnames(seu) <- paste(prefix, colnames(seu), sep = "_")
  seu_singlets_list[[i]] <- seu
}

# sanity check: no duplicates across all cells
all_cells <- unlist(lapply(seu_singlets_list, colnames))
anyDuplicated(all_cells)  # should be 0

# Skipping merge/integrate layers for v4 workflow ----
# Use the list of singlet objects instead of merged
seu_singlets_list <- lapply(
  seu_singlets_list,
  function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x)
    x <- ScaleData(x)
    x <- RunPCA(x)
    x
  }
)

anchors <- FindIntegrationAnchors(
  object.list = seu_singlets_list,
  reduction   = "rpca"  # or "cca"
)

seu_integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(seu_integrated) <- "integrated"
seu_integrated <- ScaleData(seu_integrated)
seu_integrated <- RunPCA(seu_integrated)
seu_integrated <- FindNeighbors(seu_integrated, dims = 1:20)
seu_integrated <- FindClusters(seu_integrated, resolution = 0.03)
seu_integrated <- RunUMAP(seu_integrated, dims = 1:20)

saveRDS(seu_integrated, file = "seu_integrated.rds")