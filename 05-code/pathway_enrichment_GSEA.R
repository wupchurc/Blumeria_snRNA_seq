# ---- Import Packages and Data ----
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

# ---- Combined GSEA ----
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
    fun = "gseGO",
    OrgDb = org.Rn.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
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
    
  p_gsea <- dotplot(comp_gsea, showCategory = 5, font.size = 8,
                      split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste(cell_name, "- GO (BP)"))
    
  print(p_gsea)
    
  list(gsea = comp_gsea, ranked_lists = ranked_lists)
}

# ---- Single GSEA ----
run_single_gsea <- function(ranked_vec) {
  gseGO(
    geneList = ranked_vec,
    OrgDb = org.Rn.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
}

make_single_gsea_list <- function(gsea_obj) {
  list(
    water_vs_ctrl = run_single_gsea(gsea_obj$ranked_lists$water_vs_ctrl),
    water_vs_blum = run_single_gsea(gsea_obj$ranked_lists$water_vs_blum)
  )
}

# ---- Cardiomyocytes ----
gsea_res <- analyze_celltype_gsea(cm_results, "Cardiomyocytes")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)
write.csv(
  gsea_df,
  file = "04-results/CM_GSEA_GO.csv",
  row.names = FALSE
)

# ---- Fibroblasts ----

gsea_res <- analyze_celltype_gsea(fb_results, "Fibroblasts")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)
write.csv(
  gsea_df,
  file = "04-results/FB_GSEA_GO.csv",
  row.names = FALSE
)

# ---- Macrophages ----

gsea_res <- analyze_celltype_gsea(mac_results, "Macrophages")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)
write.csv(
  gsea_df,
  file = "04-results/MAC_GSEA_GO.csv",
  row.names = FALSE
)



# ----
# Simplify the combined result
gsea_simplified <- simplify(gsea_res$gsea, cutoff = 0.7, by = "p.adjust", 
                            measure = "Wang")
gsea_simplified_df <- as.data.frame(gsea_simplified)
write.csv(
  gsea_simplified_df,
  file = "04-results/CM_GSEA_compareCluster_results_simplified.csv",
  row.names = FALSE
)

dotplot(gsea_simplified, showCategory = 5, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("CM GSEA - Top Terms") +
  theme(axis.text.y = element_text(size = 9))


# 
dotplot(gsea_res$gsea, showCategory = 5, split = ".sign", size = "Count") + 
  facet_grid(. ~ .sign) + 
  ggtitle("Pericytes - GSEA") +
  theme(axis.text.y = element_text(size = 8))


# ---- single plots ---- 
gsea_single <- make_single_gsea_list(gsea_res)
# Convert ENTREZID → SYMBOL
gsea_single$water_vs_ctrl <- setReadable(gsea_single$water_vs_ctrl, 
                                         OrgDb = org.Rn.eg.db, 
                                         keyType = "ENTREZID")
gsea_single$water_vs_blum <- setReadable(gsea_single$water_vs_blum, 
                                         OrgDb = org.Rn.eg.db, 
                                         keyType = "ENTREZID")
# emapplot
# Compute similarities for ctrl
gsea_single_ctrl_sim <- pairwise_termsim(gsea_single$water_vs_ctrl, method = "JC")
emap_ctrl <- emapplot(gsea_single_ctrl_sim, showCategory = 20, layout = "nicely")
print(emap_ctrl)
# Same for blum
gsea_single_blum_sim <- pairwise_termsim(gsea_single$water_vs_blum, method = "JC")
emap_blum <- emapplot(gsea_single_blum_sim, showCategory = 20, layout = "nicely")
print(emap_blum)

# cnetplot
cnet_ctrl <- cnetplot(gsea_single$water_vs_ctrl, 
                      showCategory = 5,
                      # categorySize = 1.2,
                      node_label = "share",     
                      layout = "fr")   
print(cnet_ctrl)
cnet_blum <- cnetplot(gsea_single$water_vs_blum, 
                      showCategory = 5,
                      # categorySize = 1.2,
                      node_label = "share",     
                      layout = "fr") 
print(cnet_blum)

#
dotplot(gsea_single$water_vs_ctrl, showCategory = 30) + ggtitle("CM water_vs_ctrl")
dotplot(gsea_single$water_vs_blum, showCategory = 10) + ggtitle("CM water_vs_blum")



