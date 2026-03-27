library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(DESeq2)
library(ggtangle)
library(ggplot2)

cm_results <- readRDS("03-analysis_scratch/DEG_Cardiomyocytes.rds")
fb_results <- readRDS("03-analysis_scratch/DEG_Fibroblasts.rds")
mac_results <- readRDS("03-analysis_scratch/DEG_Macrophages.rds")
nphl_results <- readRDS("03-analysis_scratch/DEG_Neutrophils.rds")
peri_results <- readRDS("03-analysis_scratch/DEG_Pericytes.rds")
bcells_results <- readRDS("03-analysis_scratch/DEG_B_cells.rds")
cap_results <- readRDS("03-analysis_scratch/DEG_Capillary_EC.rds")
lymph_results <- readRDS("03-analysis_scratch/DEG_Lymphatic_EC.rds")
neuro_results <- readRDS("03-analysis_scratch/DEG_Neuronal.rds")
tcell_results <- readRDS("03-analysis_scratch/DEG_T_cells.rds")
ven_results <- readRDS("03-analysis_scratch/DEG_Venous_EC.rds")

analyze_celltype_gsea <- function(results_obj, cell_name) {
  
  # Function to get ranked list
  get_ranked_list <- function(res) {
    df <- as.data.frame(res)
    df$gene <- rownames(res)
    id_map <- bitr(df$gene, fromType="SYMBOL", 
                   toType="ENTREZID", OrgDb=org.Rn.eg.db)
    df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
    df <- df[!is.na(df$ENTREZID) & !is.na(df$log2FoldChange), ]
    gene_list <- df$log2FoldChange
    names(gene_list) <- df$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)
    return(gene_list)
  }
  
  ranked_lists <- list(
    water_vs_ctrl = get_ranked_list(results_obj$water_vs_ctrl),
    water_vs_blum = get_ranked_list(results_obj$water_vs_blum)
  )
  
  gsea_df <- do.call(rbind, lapply(names(ranked_lists), function(nm) {
    data.frame(
      ENTREZID = names(ranked_lists[[nm]]),
      log2FoldChange = ranked_lists[[nm]],
      contrast_name = nm,
      stringsAsFactors = FALSE
    )
  }))
  # Set desired order here
  gsea_df$contrast_name <- factor(gsea_df$contrast_name,
                                  levels = c("water_vs_ctrl", "water_vs_blum"))
  
  comp_gsea <- compareCluster(
    ENTREZID | log2FoldChange ~ contrast_name,
    data = gsea_df,
    fun = "gseKEGG",
    organism = "rno",           # rat KEGG code
    keyType = "ncbi-geneid",    # KEGG expects this for NCBI IDs
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
  
  # Make sure Cluster is ordered consistently for plotting
  comp_gsea@compareClusterResult$Cluster <- factor(
    comp_gsea@compareClusterResult$Cluster,
    levels = c("water_vs_ctrl", "water_vs_blum")
  )
  
  list(gsea = comp_gsea, ranked_lists = ranked_lists)
}


#---- Cardiomyocytes ----
gsea_res <- analyze_celltype_gsea(cm_results, "Cardiomyocytes")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/CM_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Cardiomyocytes - KEGG") +
  theme(axis.text.y = element_text(size = 9))

#----
gsea_res <- analyze_celltype_gsea(mac_results, "Macrophages")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/MAC_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Macrophages - KEGG") +
  theme(axis.text.y = element_text(size = 9))

# 
gsea_res <- analyze_celltype_gsea(fb_results, "Fibroblasts")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/FB_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Fibroblasts - KEGG") +
  theme(axis.text.y = element_text(size = 9))


#---- Neutrophils ----
gsea_res <- analyze_celltype_gsea(nphl_results, "Neutrophils")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/NPHL_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Neutrophils - KEGG") +
  theme(axis.text.y = element_text(size = 9))

#---- B cells ----
gsea_res <- analyze_celltype_gsea(bcells_results, "B cells")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/BCELLS_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("B cells - KEGG") +
  theme(axis.text.y = element_text(size = 9))
#---- Pericytes ----
gsea_res <- analyze_celltype_gsea(peri_results, "Pericytes")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/PERI_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Pericytes - KEGG") +
  theme(axis.text.y = element_text(size = 9))

# ---- Capillary EC ----
gsea_res <- analyze_celltype_gsea(cap_results, "Capillary EC")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/CAPILLARYEC_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Capillary EC - KEGG") +
  theme(axis.text.y = element_text(size = 9))


#---- Lymphatic EC ----
gsea_res <- analyze_celltype_gsea(lymph_results, "Lymphatic EC")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/LYMPHATICEC_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Lymphatic EC - KEGG") +
  theme(axis.text.y = element_text(size = 9))
# ---- Neuronal ----
gsea_res <- analyze_celltype_gsea(neuro_results, "Neuronal")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/NEURONAL_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Neuronal - KEGG") +
  theme(axis.text.y = element_text(size = 9))
# ---- T cells ----
gsea_res <- analyze_celltype_gsea(tcell_results, "T cells")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/TCELL_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("T cells - KEGG") +
  theme(axis.text.y = element_text(size = 9))

# ---- Venous EC ----
gsea_res <- analyze_celltype_gsea(ven_results, "Venous EC")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)

write.csv(
  gsea_df,
  file = "04-results/VENOUSEC_GSEA_KEGG_results.csv",
  row.names = FALSE
)

dotplot(gsea_res$gsea, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Venous EC - KEGG") +
  theme(axis.text.y = element_text(size = 9))
