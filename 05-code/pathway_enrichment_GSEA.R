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


analyze_celltype_gsea <- function(results_obj, cell_name) {  # No lfc_threshold needed
  
  # Function to get ranked list
  get_ranked_list <- function(res) {
    df <- as.data.frame(res)
    df$gene <- rownames(res)
    id_map <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)
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
  
  # GSEA comparison with formula syntax for compareCluster
  comp_gsea <- compareCluster(
    ENTREZID | log2FoldChange ~ contrast_name,
    data = do.call(rbind, lapply(names(ranked_lists), function(nm) {
      data.frame(ENTREZID = names(ranked_lists[[nm]]), 
                 log2FoldChange = ranked_lists[[nm]], 
                 contrast_name = nm)
    })),
    fun = "gseGO",
    OrgDb = org.Rn.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    nPermSimple = 10000,  # Permutations for GSEA
    verbose = FALSE
  )
  
  # Plot
  p_gsea <- dotplot(comp_gsea, showCategory = 5, font.size = 8, split = ".sign") + 
    facet_grid(. ~ .sign) +
    ggtitle(paste(cell_name, "- GSEA (BP)"))
  
  print(p_gsea)
  
  return(list(gsea = comp_gsea, ranked_lists = ranked_lists))
}


cm_gsea <- analyze_celltype_gsea(cm_results, "Cardiomyocytes")
fb_gsea <- analyze_celltype_gsea(fb_results, "Fibroblasts")
mac_gsea <- analyze_celltype_gsea(mac_results, "Macrophages")
nphl_gsea <- analyze_celltype_gsea(nphl_results, "Neutrophils")
peri_gsea <- analyze_celltype_gsea(peri_results, "Pericytes")