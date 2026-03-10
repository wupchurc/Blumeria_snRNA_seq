library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(DESeq2)
library(ggtangle)
library(ggplot2)

cm_results <- readRDS("03-analysis_scratch/DEG_Cardiomyocytes.rds")
mac_results <- readRDS("03-analysis_scratch/DEG_Macrophages.rds")

#---- Function to process one DESeqResults → up/down ENTREZ lists ----
get_sig_genes <- function(res, lfc_threshold = 0.5, padj_threshold = 0.05, 
                          comp_name = NULL, cell_name = NULL, write_files = FALSE) {
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  
  id_map <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)
  df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
  
  up   <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange >  lfc_threshold]
  down <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange < -lfc_threshold]
  
  # Write files if requested
  if (write_files && !is.null(comp_name) && !is.null(cell_name)) {
    
    up_genes <- up[!is.na(up)]
    down_genes <- down[!is.na(down)]
    background <- as.character(df$ENTREZID)
    
    writeLines(up_genes, sprintf("03-analysis_scratch/%s_%s_up_genes.txt", cell_name, comp_name))
    writeLines(down_genes, sprintf("03-analysis_scratch/%s_%s_down_genes.txt", cell_name, comp_name))
    writeLines(background, sprintf("03-analysis_scratch/%s_%s_background.txt", cell_name, comp_name))
  }
  
  
  list(up = up, down = down, all = df$ENTREZID)
}
# ---- Print up/down regulated DEGs and background to txt file (can be used w/ ShinyGO) ----

get_sig_genes(res = cm_results$water_vs_ctrl, lfc_threshold = 0.5, padj_threshold = 0.05, 
              comp_name = "water_vs_ctrl", cell_name = "cm", write_files = TRUE)
get_sig_genes(res = cm_results$water_vs_blum, lfc_threshold = 0.5, padj_threshold = 0.05, 
              comp_name = "water_vs_blum", cell_name = "cm", write_files = TRUE)


# ---- Refactored Code for ORA ----

analyze_celltype <- function(results_obj, cell_name, lfc_threshold = 0.5, padj_threshold = 0.05) {
  
  sig_lists <- list(
    water_vs_ctrl = get_sig_genes(results_obj$water_vs_ctrl, lfc_threshold, padj_threshold),
    water_vs_blum = get_sig_genes(results_obj$water_vs_blum, lfc_threshold, padj_threshold)
  )
  
  #Shared background
  shared_background <- unique(c(
    sig_lists$water_vs_ctrl$all,
    sig_lists$water_vs_blum$all
  ))
  
  # Compare upregulated
  comp_up <- compareCluster(
    geneClusters = lapply(sig_lists, '[[', "up"),
    fun = "enrichGO",
    # universe = shared_background,
    universe = shared_background,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 2,
    readable = TRUE
  )
  
  # Compare downregulated
  comp_down <- compareCluster(
    geneClusters = lapply(sig_lists, '[[', "down"),
    fun = "enrichGO",
    universe = shared_background,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 2,
    readable = TRUE
  )
  
  # Plots
  p_up <- dotplot(comp_up, showCategory = 10, font.size = 8, size = "Count") + 
    ggtitle(paste(cell_name, "- Upregulated (p<0.25)")) 
  # +
    # scale_size_continuous(name = "Gene\nCount", breaks = c(5,10,15,20,25,30,35,40,45,50),
                          # limits = c(0,50), range = c(2,10))
  
  p_down <- dotplot(comp_down, showCategory = 10, font.size = 8, size = "Count") +
    ggtitle(paste(cell_name, "- Downregulated (p<0.25)")) 
  # +
    # scale_size_continuous(name = "Gene\nCount", breaks = c(5,10,15,20,25,30,40,45,50),
                          # limits = c(0,50), range = c(2,10))
  
  print(p_up)
  print(p_down)
  
  # Return results
  list (up = comp_up, down = comp_down, sig_lists = sig_lists, background = shared_background)
  
}



cm_analysis <- analyze_celltype(cm_results, "Cardiomyocytes")



# ---- GSEA ----
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
mac_gsea <- analyze_celltype_gsea(mac_results, "Macrophages")
