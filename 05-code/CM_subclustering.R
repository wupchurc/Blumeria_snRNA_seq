
#If you don't already have the cm specific data: 
#CM_data  <- subset(seu_integrated, idents = c("Cardiomyocytes"))

#If you do: 
CM_data <- readRDS("/projects/standard/szprisco/shared/murp1832_projects/2026/CM_subclustering/Blumeria_scRNA_seq/03-analysis_scratch/CM_data.rds")

CM_data <- FindNeighbors(CM_data, reduction = "harmony", dims = 1:20)
CM_data <- FindClusters(CM_data, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
CM_data <- RunUMAP(CM_data, reduction = "harmony", dims = 1:20)

CM_data$condition <- factor(
  CM_data$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)
DimPlot(
  CM_data,
  reduction = "umap",
  group.by = "condition",
  split.by = "condition",
  label = FALSE
)

DimPlot(
  CM_data,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  split.by = "condition",
  label = TRUE
)

umap_1 <- DimPlot(
  CM_data,
  reduction = "umap",
  group.by = "condition",
  split.by = "condition",
  label = FALSE
)

umap_2 <- DimPlot(
  CM_data,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  split.by = "condition",
  label = TRUE
)

ggsave(filename = "04-results/UMAP_1.png", plot = umap_1)

ggsave(filename = "04-results/UMAP_2.png", plot = umap_2)

# use RNA counts/normalized data for DE
DefaultAssay(CM_data) <- "RNA"
CM_data <- JoinLayers(CM_data, assay = "RNA")

# set identities to clusters (from integrated clustering)
Idents(CM_data) <- "seurat_clusters"

# ---- conserved markers per cluster across conditions ----
# Run FindConservedMarkers on clusters to help identify
markers_cluster0  <- FindConservedMarkers(CM_data, ident.1 = 0,  grouping.var = "condition")
markers_cluster1  <- FindConservedMarkers(CM_data, ident.1 = 1,  grouping.var = "condition")
markers_cluster2  <- FindConservedMarkers(CM_data, ident.1 = 2,  grouping.var = "condition")
markers_cluster3  <- FindConservedMarkers(CM_data, ident.1 = 3,  grouping.var = "condition")
markers_cluster4  <- FindConservedMarkers(CM_data, ident.1 = 4,  grouping.var = "condition")

# Store renamed Idents in metadata column
CM_data$cell_type <- Idents(CM_data)

# Set the desired order of cell types
celltype_order <- c(
  "0",
  "1",
  "2",
  "3",
  "4"
)

# Define colors matching celltype_order 
cell_cols <- c(
  "0" = "#BA1C30",
  "1"    = "#5FA641", 
  "2"      = "#702C8C",
  "3"      = "#CC79A7",
  "4"   = "#D55E00"
)

# Make cell_type a factor with that order
CM_data$cell_type <- factor(
  CM_data$cell_type,
  levels = celltype_order
)

# --- plotting (still using integrated UMAP) ----
DefaultAssay(CM_data) <- "integrated"  # for plotting/clustering context

# Reorder condition for plots
CM_data$condition <- factor(
  CM_data$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)

# Plot UMAPs with identified clusters
print(DimPlot(CM_data, reduction = 'umap', 
              label = TRUE, label.size = 4, repel = TRUE, cols = cell_cols) + NoLegend())

# ---- Relative quantification ----

# After setting celltype_order, reverse it for plotting consistency
plot_order <- rev(celltype_order)  # Reverse for ggplot/barplot
# Set Idents and cell_type with plot order
Idents(CM_data) <- factor(CM_data$cell_type, levels = plot_order)
CM_data$cell_type <- factor(CM_data$cell_type, levels = plot_order)
# Now cell_types and plots follow your desired order
cell_types <- levels(Idents(CM_data))

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
  CM_data,
  reduction  = "umap",
  split.by   = "condition",
  label      = TRUE,
  label.size = 4,
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
  CellType  = Idents(CM_data),
  Condition = CM_data$condition
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

# ---- sample-level Principal Component Analysis ----

avg <- AverageExpression(
  CM_data, 
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

df$condition <- CM_data$condition[match(df$sample,
                                        CM_data$orig.ident)]
condition_levels <- levels(CM_data$condition)

p_pca <- ggplot(df, aes(PC1, PC2, color = condition, fill = condition)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, alpha = 0.3, geom = "polygon") +
  labs(x = paste0("PC1 (", pc_var[1], "%)"),
       y = paste0("PC2 (", pc_var[2], "%)")) + 
  theme_classic() +
  theme(legend.title = element_blank()) # Removes legend title

print(p_pca)

CM_obj <- readRDS("/projects/standard/szprisco/shared/murp1832_projects/2026/CM_subclustering/Blumeria_scRNA_seq/03-analysis_scratch/CM_for_DGE.rds")

#Consolidated Pseudobulking and DESeq2 Function ----
run_pseudobulk_deg <- function(seu_obj, cell_type, min_counts = 10, alpha = 0.05,
                               plot_title = NULL, save_results = FALSE) {
  
  # Set assay to RNA for raw counts
  DefaultAssay(seu_obj) <- "RNA"
  
  # Subset to specified cell type
  cells <- WhichCells(seu_obj, idents = cell_type)
  if (length(cells) == 0) {
    stop(paste("No cells found for cell type:", cell_type))
  }
  
  seu_subset <- subset(seu_obj, cells = cells)
  cat("Analyzing", length(cells), "cells for", cell_type, "\n") # look at, maybe remove
  
  # Pseudobulk Matrix by sample (orig.ident)
  cts <- AggregateExpression(seu_subset, assays = "RNA", group.by = "orig.ident",
                             slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  
  # Prepare Metadata
  colData <- data.frame(samples = colnames(cts)) %>%
    mutate(
      condition = case_when(
        str_detect(samples, "Blumeria") ~ "MCT-Blumeria",
        str_detect(samples, "Water") ~ "MCT-Water",
        str_detect(samples, "Control") ~ "Control"
      )
    ) %>%
    column_to_rownames(var = "samples")
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = colData,
    design = ~ condition)
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= min_counts
  dds <- dds[keep, ]
  cat("Filtered to", nrow(dds), "genes (min total counts =", min_counts, ")\n")
  
  # Set reference level 
  dds$condition <- relevel(dds$condition, ref = "Control")
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract contrasts
  res_water_vs_ctrl <- results(dds, contrast = c("condition", "MCT-Water", "Control"),
                               alpha = alpha)
  res_blum_vs_ctrl <- results(dds, contrast = c("condition", "MCT-Blumeria", "Control"),
                              alpha = alpha)
  res_blum_vs_water <- results(dds, contrast = c("condition", "MCT-Blumeria", "MCT-Water"), 
                               alpha = alpha)
  res_water_vs_blum <- results(dds, contrast = c("condition", "MCT-Water", "MCT-Blumeria"),
                               alpha = alpha)
  
  # Print summaries
  cat("\n=== DESeq2 Results Summary for", cell_type, "===\n")
  cat("Water vs Control:\n"); summary(res_water_vs_ctrl)
  cat("\nBlumeria vs Control:\n"); summary(res_blum_vs_ctrl)
  cat("\nBlumeria vs Water:\n"); summary(res_blum_vs_water)
  cat("\nWater vs Blumeria:\n"); summary(res_water_vs_blum)
  
  # Normalize data for PCA
  # rld <- vst(dds, blind = FALSE)
  
  # Generate plots
  # title <- ifelse(is.null(plot_title), cell_type, plot_title)
  
  # p_pca <- plotPCA(rld, intgroup = "condition") +
  # stat_ellipse(level = 0.95, alpha = 0.3, geom = "polygon") +
  # ggtitle(paste("Pseudoulk PCA:", title))
  
  # p_ma_water <- plotMA(res_water_vs_ctrl) + ggtitle("MA: Water vs Control")
  # p_ma_water_blum <- plotMA(res_water_vs_blum) + ggtitle("MA: Water vs Blumeria")
  
  # Display plots
  # print(p_pca)
  # print(p_ma_water)
  # print(p_ma_water_blum)
  
  # Save results if requested
  if (save_results) {
    results_list <- list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      water_vs_blum = res_water_vs_blum
      # dds = dds
      # rld = rld
    )
    
    filename <- paste0("DEG_", gsub("[^A-Za-z0-9]", "_", cell_type), ".rds")
    saveRDS(results_list, file = paste0("/projects/standard/szprisco/shared/murp1832_projects/2026/CM_subclustering/Blumeria_scRNA_seq/03-analysis_scratch/", filename))
    cat("Results saved to:", filename, "\n")
    
  }
  
  # Return results
  return(
    results <- list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      water_vs_blum = res_water_vs_blum,
      dds = dds
      # rld = rld
    )
  )
}

# ---- Bar Plots of up and down regulated gene counts ----
# Run DSEq2 on all cell types
deg_results <- list()
cell_types <- levels(CM_obj)  
for (cell_type in cell_types) {
  cat("\\\\n=== Processing", cell_type, "===\\\\n")
  deg_results[[cell_type]] <- run_pseudobulk_deg(CM_obj, cell_type, 
                                                 alpha = 0.1, save_results = TRUE)
}








