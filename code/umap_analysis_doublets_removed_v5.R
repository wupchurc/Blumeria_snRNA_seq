# load libraries ----
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(presto)
library(patchwork)
library(scCustomize)

# Load control dataset from 10X CellRanger ----
data_dirs <- c("data/Analysis/cellranger_Control_SP_11_GEX_FL-Z0041/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_Control_SP_12_GEX_FL-Z0044/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_Control_SP_14_GEX_FL-Z0047/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_Control_SP_15_GEX_FL-Z0054/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_1_GEX_FL-Z0043/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_5_GEX_FL-Z0046/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_6_GEX_FL-Z0053/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Water_9_GEX_FL-Z0056/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_1_GEX_FL-Z0042/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_3_GEX_FL-Z0045/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_4_GEX_FL-Z0048/raw_feature_bc_matrix/",
               "data/Analysis/cellranger_MCT_Blumeria_7_GEX_FL-Z0055/raw_feature_bc_matrix/")
sample_ids <- c("Control11", "Control12", "Control14", "Control15",
                "Water1", "Water5", "Water6", "Water9",
                "Blumeria1", "Blumeria3", "Blumeria4", "Blumeria7")

# Create seurat objects ----
obj_list <- vector("list", length(data_dirs))

for (i in seq_along(data_dirs)) {
  mtx <- Read10X(data.dir = data_dirs[i])
  seu <- CreateSeuratObject(
    counts       = mtx,
    min.cells    = 3,
    min.features = 200,
    project      = sample_ids[i]
  )
  seu$sample    <- sample_ids[i]
  seu$condition <- dplyr::case_when(
    grepl("^Control",  sample_ids[i]) ~ "Control",
    grepl("^Water",    sample_ids[i]) ~ "MCT-Water",
    grepl("^Blumeria", sample_ids[i]) ~ "MCT-Blumeria"
  )
  obj_list[[i]] <- seu
}

# Doublets ----
# Doublet rate estimation
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
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
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
# ----

# saveRDS(doubletFinder_obj_list, file = "doubletFinder_obj.rds")
doubletFinder_obj_list <- readRDS("doubletFinder_obj.rds")

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

# 2.1 Merge all singlets
seu_merged <- merge(
  x = seu_singlets_list[[1]],
  y = seu_singlets_list[-1],
  #add.cell.ids = sample_ids  # you already used prefixing, this is optional
)

DefaultAssay(seu_merged) <- "RNA"
# create a 'counts' layer (base) if not already present
Layers(seu_merged)

# START OVER from after doublet removal / merge (where you had counts.* layers)
# Skip JoinLayers entirely!

DefaultAssay(seu_merged) <- "RNA"

# Clean normalization WITHOUT joining
seu_merged <- NormalizeData(seu_merged, verbose = FALSE)  # normalizes each layer

# Find variable features on shared genes only
count_layers <- grep("^counts\\.", Layers(seu_merged), value = TRUE)
shared_features <- Reduce(intersect, lapply(count_layers, function(ly) rownames(LayerData(seu_merged, layer = ly))))

# seu_merged <- FindVariableFeatures(seu_merged, 
                                   # features = shared_features, 
                                   # nfeatures = 2000)

# Scale/PCA (layer-aware)
seu_merged <- ScaleData(seu_merged, verbose = FALSE)
seu_merged <- RunPCA(seu_merged, npcs = 30, verbose = FALSE)

# NOW IntegrateLayers with layers= argument
sample_ids <- gsub("^counts\\.", "", count_layers)
# NOW IntegrateLayers - use metadata (no layers= needed)
seu_merged <- IntegrateLayers(
  object        = seu_merged,
  assay         = "RNA",
  method        = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  layers = count_layers
)


saveRDS(seu_merged, file = "seu_merged_integrated.rds")
# seu_merged <- readRDS("seu_merged_integrated.rds")

