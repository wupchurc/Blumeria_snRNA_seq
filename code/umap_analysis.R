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

# Merge 
merged_seurat <- merge(obj_list[[1]], y = obj_list[-1], 
                         add.cell.ids = sample_ids)

# QC & Filtering ----

print(FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
        geom_smooth(method = 'lm'))

# calculate % mt
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# filter
merged_seurat <- subset(merged_seurat,
                        subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Perform standard workflow steps to first visualize batch effects ----
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(object = merged_seurat)
merged_seurat <- ScaleData(object = merged_seurat)
merged_seurat <- RunPCA(object = merged_seurat)
print(ElbowPlot(merged_seurat))
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:15)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5,
                              cluster.name = "unintegrated_clusters")
merged_seurat <- RunUMAP(merged_seurat, dims = 1:15, 
                         reduction.name = "umap.unintegrated")

print(DimPlot(merged_seurat, reduction = "umap.unintegrated", 
        group.by = c("condition", "unintegrated_clusters")))
# Integrate ----
merged_seurat <- IntegrateLayers(object = merged_seurat, method = CCAIntegration,
                                 orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)

#save/load integrated object
#saveRDS(merged_seurat, file = "integrated_object.rds")
merged_seurat <- readRDS("integrated_object.rds")

# re-join layers after integration
merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
#
merged_seurat <- FindNeighbors(merged_seurat, reduction = "integrated.cca", dims = 1:15)
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.025,0.05,0.075,0.1))
print(DimPlot(merged_seurat, reduction = 'pca', 
              group.by='RNA_snn_res.0.025',label=TRUE))
Idents(merged_seurat) <- 'RNA_snn_res.0.025'
merged_seurat <- RunUMAP(merged_seurat, dims = 1:15, reduction = "integrated.cca")


print(DimPlot(merged_seurat, reduction = 'umap', label = TRUE))
print(DimPlot(merged_seurat, reduction = "umap", split.by  = "condition",
              group.by = 'condition', ncol = 3) + NoLegend() + ggtitle(NULL))

# Load integrated object for cluster identification ----
#saveRDS(merged_seurat, file = "merged_seurat_obj.rds")
merged_seurat <- readRDS("merged_seurat_obj.rds")

# findConserved markers ----
markers_cluster0 <- FindConservedMarkers(merged_seurat, ident.1 = 0,
                                         grouping.var = 'condition')
markers_cluster1 <- FindConservedMarkers(merged_seurat, ident.1 = 1,
                                         grouping.var = 'condition')
markers_cluster2 <- FindConservedMarkers(merged_seurat, ident.1 = 2,
                                         grouping.var = 'condition')
markers_cluster3 <- FindConservedMarkers(merged_seurat, ident.1 = 3,
                                         grouping.var = 'condition')
markers_cluster4 <- FindConservedMarkers(merged_seurat, ident.1 = 4,
                                         grouping.var = 'condition')
markers_cluster5 <- FindConservedMarkers(merged_seurat, ident.1 = 5,
                                         grouping.var = 'condition')
markers_cluster6 <- FindConservedMarkers(merged_seurat, ident.1 = 6,
                                         grouping.var = 'condition')
markers_cluster7 <- FindConservedMarkers(merged_seurat, ident.1 = 7,
                                         grouping.var = 'condition')
markers_cluster8 <- FindConservedMarkers(merged_seurat, ident.1 = 8,
                                         grouping.var = 'condition')
markers_cluster9 <- FindConservedMarkers(merged_seurat, ident.1 = 9,
                                         grouping.var = 'condition')
markers_cluster10 <- FindConservedMarkers(merged_seurat, ident.1 = 10,
                                          grouping.var = 'condition')

print(FeaturePlot(merged_seurat, reduction='umap', features = c("Son", "Tnrc6b", "Pdgfra"),
                  min.cutoff = 'q10'))
print(FeaturePlot(merged_seurat, reduction='umap', features = c("Son", "Tnrc6b", "Ldb3"),
                  min.cutoff = 'q10'))
print(FeaturePlot(merged_seurat, reduction = 'umap', features = "Cd79a"))

print(FeaturePlot(merged_seurat, reduction = 'umap', features = c('Lyve1','Gja5')))
# Renaming clusters ----
# rename cluster 0
print(Idents(merged_seurat))
merged_seurat <- RenameIdents(merged_seurat, '0' = 'Fibroblasts')
# rename cluster 1
merged_seurat <- RenameIdents(merged_seurat, '1' = 'Capillary EC')
# rename cluster 2
merged_seurat <- RenameIdents(merged_seurat, '2' = 'Macrophages')
# rename cluster 3
merged_seurat <- RenameIdents(merged_seurat, '3' = 'Cardiomyocytes')
# rename cluster 4
merged_seurat <- RenameIdents(merged_seurat, '4' = 'Doublets*')
# rename cluster 5
merged_seurat <- RenameIdents(merged_seurat, '5' = 'Mural cells')
# rename cluster 6
merged_seurat <- RenameIdents(merged_seurat, '6' = 'T cells')
# rename cluster 7
merged_seurat <- RenameIdents(merged_seurat, '7' = 'B cells')
# rename cluster 8 
merged_seurat <- RenameIdents(merged_seurat, '8' = 'Lymphatic EC')
# rename cluster 9
merged_seurat <- RenameIdents(merged_seurat, '9' = 'Cardiac neurons')
# rename cluster 10
merged_seurat <- RenameIdents(merged_seurat, '10' = 'Neutrophils')

merged_seurat$condition <- factor(merged_seurat$condition,
                                  levels = c("Control", "MCT-Water", "MCT-Blumeria"))
# Replot UMAPs with identified clusters
print(DimPlot(merged_seurat, reduction = 'umap',split.by = 'condition', label = TRUE, 
              label.size = 3, repel      = TRUE)) + NoLegend()
print(DimPlot(merged_seurat, reduction = 'umap',split.by = 'condition', label = FALSE)) 

# Relative quantification ----

# shared colors (same order as Idents)
cell_types <- levels(Idents(merged_seurat))
#cell_cols  <- Seurat::DiscretePalette(length(cell_types))
cell_cols <- DiscretePalette_scCustomize(
  num_colors = length(cell_types),
  palette    = "glasbey"   # or "polychrome", "alphabet2", "varibow"
)
names(cell_cols) <- cell_types


origin_x <- -15
origin_y <- -15
len      <- 5   # length of mini-axes in UMAP units

mini_axes <- data.frame(
  x    = c(origin_x, origin_x),
  y    = c(origin_y, origin_y),
  xend = c(origin_x + len, origin_x),
  yend = c(origin_y,        origin_y + len)
)
# Use same palette in UMAP
p_umap <- DimPlot(
  merged_seurat,
  reduction  = "umap",
  split.by   = "condition",
  label      = TRUE,
  label.size = 3,
  cols       = cell_cols          # <— key line
) +
  NoGrid() + 
  NoLegend() +
  theme(
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    axis.text   = element_blank(),
    axis.title  = element_blank(),
  ) +
  geom_segment(
    data = mini_axes,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(4, "pt")),
    linewidth = 0.4
  ) +
  annotate("text",
           x = origin_x + len, y = origin_y,
           label = "UMAP1", vjust = 1.5, size = 3) +
  annotate("text",
           x = origin_x - 0.25, y = origin_y + len,
           label = "UMAP2", hjust = 0.5, angle = 90, size = 3)

# counts per cell type per condition
tab_counts <- table(
  CellType  = Idents(merged_seurat),
  Condition = merged_seurat$condition
)

# relative proportions within each condition (columns sum to 1)
tab_prop <- prop.table(tab_counts, margin = 2)
df_prop <- as.data.frame(tab_prop) %>%
  rename(Proportion = Freq)

p_bar <- ggplot(df_prop,
                aes(x = 1, y = Proportion, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 1) +
  coord_flip() +
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  labs(title = "Relative Abundance (%)") +   # ← add title here
  xlab(NULL) +   # now appears above the bars
  ylab(NULL) +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 10),  # center title
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),  # keep digits hidden
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_blank(),   
    axis.ticks.y     = element_blank()
  )

# plot combined graphs
p_umap / p_bar + plot_layout(heights = c(12, 1))

