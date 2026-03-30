# load libraries 
library(Seurat)
library(DESeq2)
library(tidyverse)
library(patchwork)

library(presto)
library(scCustomize)
library(SummarizedExperiment)
library(RColorBrewer)
library(circlize)
library(tidyr)

# load processed seurat object from rds file 
seu_integrated <- readRDS("03-analysis_scratch/seu_integrated.rds")

# use RNA counts/normalized data for DE
DefaultAssay(seu_integrated) <- "RNA"
seu_integrated <- JoinLayers(seu_integrated, assay = "RNA")

# set identities to clusters (from integrated clustering)
Idents(seu_integrated) <- "seurat_clusters"

# ---- conserved markers per cluster across conditions ----
# Run FindConservedMarkers on clusters to help identify
markers_cluster0  <- FindConservedMarkers(seu_integrated, ident.1 = 0,  grouping.var = "condition")
markers_cluster1  <- FindConservedMarkers(seu_integrated, ident.1 = 1,  grouping.var = "condition")
markers_cluster2  <- FindConservedMarkers(seu_integrated, ident.1 = 2,  grouping.var = "condition")
markers_cluster3  <- FindConservedMarkers(seu_integrated, ident.1 = 3,  grouping.var = "condition")
markers_cluster4  <- FindConservedMarkers(seu_integrated, ident.1 = 4,  grouping.var = "condition")
markers_cluster5  <- FindConservedMarkers(seu_integrated, ident.1 = 5,  grouping.var = "condition")
markers_cluster6  <- FindConservedMarkers(seu_integrated, ident.1 = 6,  grouping.var = "condition")
markers_cluster7  <- FindConservedMarkers(seu_integrated, ident.1 = 7,  grouping.var = "condition")
markers_cluster8  <- FindConservedMarkers(seu_integrated, ident.1 = 8,  grouping.var = "condition")
markers_cluster9  <- FindConservedMarkers(seu_integrated, ident.1 = 9,  grouping.var = "condition")
markers_cluster10 <- FindConservedMarkers(seu_integrated, ident.1 = 10, grouping.var = "condition")

# Renaming clusters 
print(Idents(seu_integrated))
seu_integrated <- RenameIdents(seu_integrated, '0' = 'Fibroblasts')
seu_integrated <- RenameIdents(seu_integrated, '1' = 'Capillary EC')
seu_integrated <- RenameIdents(seu_integrated, '2' = 'Macrophages')
seu_integrated <- RenameIdents(seu_integrated, '3' = 'Cardiomyocytes')
seu_integrated <- RenameIdents(seu_integrated, '4' = 'Pericytes')
seu_integrated <- RenameIdents(seu_integrated, '5' = 'T cells')
seu_integrated <- RenameIdents(seu_integrated, '6' = 'Venous EC')
seu_integrated <- RenameIdents(seu_integrated, '7' = 'B cells')
seu_integrated <- RenameIdents(seu_integrated, '8' = 'Lymphatic EC')
seu_integrated <- RenameIdents(seu_integrated, '9' = 'Neutrophils')
seu_integrated <- RenameIdents(seu_integrated, '10' = 'Neuronal')

# Store renamed Idents in metadata column
seu_integrated$cell_type <- Idents(seu_integrated)

# Set the desired order of cell types
celltype_order <- c(
  "Cardiomyocytes",
  "Fibroblasts",
  "Macrophages",
  "T cells",
  "B cells",
  "Neutrophils",
  "Pericytes",
  "Capillary EC",
  "Venous EC",
  "Lymphatic EC",
  "Neuronal"
)

# Define colors matching celltype_order 
cell_cols <- c(
  "Cardiomyocytes" = "#BA1C30",
  "Fibroblasts"    = "#5FA641", 
  "Macrophages"      = "#702C8C",
  "T cells"   = "#999999",
  "B cells"      = "#CC79A7",
  "Neutrophils"   = "#D55E00",
  "Pericytes"    = "#0072B2",
  "Capillary EC"    = "#F0E442",
  "Venous EC"        = "#009E73",
  "Lymphatic EC"        = "#56B4E9",
  "Neuronal"       = "#E69F00"
)

# Make cell_type a factor with that order
seu_integrated$cell_type <- factor(
  seu_integrated$cell_type,
  levels = celltype_order
)

# ---- Heatmap of markers used to validate cell type identities ----
all_markers <- c(
  "Tnnt2","Myh6","Ryr2",
  "Pdgfra", "Col1a1", "Dcn", 
  "Adgre1","Csf1r",
  "Scn7a", "Cdh19", "Chl1"
)

p_heat <- DoHeatmap(subset(seu_integrated, downsample = 1000),
          features = all_markers,
          size = 3,
          hjust = 0.5,
          vjust = - 0.5,
          angle = 0, 
          group.bar.height = 0.06,
          group.by = "cell_type",
          group.colors = cell_cols) + 
  guides(color = "none") + 
  labs(fill = "Z-score")

print(p_heat)

ggsave(
  filename = "04-results/celltype_markers_heatmap.png",  # or .png, .tiff
  plot     = p_heat,
  width    = 8,   # adjust as needed
  height   = 6,
  dpi      = 300
)

# --- plotting (still using integrated UMAP) ----
DefaultAssay(seu_integrated) <- "integrated"  # for plotting/clustering context

# Reorder condition for plots
seu_integrated$condition <- factor(
  seu_integrated$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)

# Plot UMAPs with identified clusters
print(DimPlot(seu_integrated, reduction = 'umap', 
              label = TRUE, label.size = 4, repel = TRUE, cols = cell_cols) + NoLegend())

# ---- Relative quantification ----

# After setting celltype_order, reverse it for plotting consistency
plot_order <- rev(celltype_order)  # Reverse for ggplot/barplot
# Set Idents and cell_type with plot order
Idents(seu_integrated) <- factor(seu_integrated$cell_type, levels = plot_order)
seu_integrated$cell_type <- factor(seu_integrated$cell_type, levels = plot_order)
# Now cell_types and plots follow your desired order
cell_types <- levels(Idents(seu_integrated))

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
  seu_integrated,
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
  CellType  = Idents(seu_integrated),
  Condition = seu_integrated$condition
)

# relative proportions within each condition (columns sum to 1)
tab_prop <- prop.table(tab_counts, margin = 2)
df_prop <- as.data.frame(tab_prop) %>%
  dplyr::rename(Proportion = Freq)

p_bar <- ggplot(df_prop,
                aes(x = 1, y = Proportion, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 1) +
  coord_flip() +
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  xlab(NULL) +   # now appears above the bars
  ylab(NULL) +
  theme_minimal() +
  labs(title = "Relative Abundance (%)") +
  theme(
    plot.title = element_text(hjust = 0.5),
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
p_combined <- p_umap / p_bar + plot_layout(heights = c(12, 1))

print(p_combined)

ggsave(
  filename = "04-results/umap_rel_abundance.png",
  plot     = p_combined,
  width    = 15,    # adjust as needed
  height   = 7.5,
  units = "in",
  dpi      = 300
)

# ---- sample-level Principal Component Analysis ----

avg <- AverageExpression(
  seu_integrated, 
  assays = "RNA",
  layer = "data",
  group.by = "orig.ident"
)$RNA

# Remove genes with zero variance across samples
gene_var <- apply(avg, 1, var)
avg_filtered <- avg[gene_var > 0, ]

mat <- t(avg_filtered)
pca <- prcomp(mat, scale. = TRUE)

# Get variance explained percentages
pc_var <- round(100 * summary(pca)$importance[2, 1:2], 1) #PC1, PC2 %

df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample = rownames(pca$x))

df$condition <- seu_integrated$condition[match(df$sample,
                                               seu_integrated$orig.ident)]
condition_levels <- levels(seu_integrated$condition)

p_pca <- ggplot(df, aes(PC1, PC2, color = condition, fill = condition)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, alpha = 0.3, geom = "polygon") +
  labs(x = paste0("PC1 (", pc_var[1], "%)"),
       y = paste0("PC2 (", pc_var[2], "%)")) + 
  theme_classic() +
  theme(legend.title = element_blank()) # Removes legend title

print(p_pca)

ggsave(
  filename = "04-results/sample_level_pca.png",
  plot     = p_pca,
  width    = 8,    # adjust as needed
  height   = 6,
  dpi      = 300
)

# ---- Save

saveRDS(seu_integrated, file = "03-analysis_scratch/seu_for_DGE.rds")
