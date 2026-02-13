library(clusterProfiler)
library(org.Rn.eg.db)
library(DESeq2)

cm_results <- readRDS("DEG_Cardiomyocytes.rds")  # DESeqResults

# Function to process one DESeqResults → up/down ENTREZ lists
get_sig_genes <- function(res, lfc_threshold = 1, padj_threshold = 0.05) {
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  
  id_map <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)
  df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
  
  up   <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange >  lfc_threshold]
  down <- df$ENTREZID[df$padj < padj_threshold & df$log2FoldChange < -lfc_threshold]
  
  list(up = up, down = down, all = df$ENTREZID)
}

# Extract sig lists for ALL contrasts
sig_lists <- list(
  water_vs_ctrl = get_sig_genes(cm_results$water_vs_ctrl),
  blum_vs_ctrl  = get_sig_genes(cm_results$blum_vs_ctrl),
  blum_vs_water = get_sig_genes(cm_results$blum_vs_water)
)

# Individual ORA results (up/down per contrast)
ego_up_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$up, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$up, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$up, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

ego_down_all <- list(
  water_vs_ctrl = enrichGO(gene=sig_lists$water_vs_ctrl$down, 
                           universe=sig_lists$water_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_ctrl  = enrichGO(gene=sig_lists$blum_vs_ctrl$down, 
                           universe=sig_lists$blum_vs_ctrl$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05),
  blum_vs_water = enrichGO(gene=sig_lists$blum_vs_water$down, 
                           universe=sig_lists$blum_vs_water$all,
                           OrgDb=org.Rn.eg.db, keyType="ENTREZID", ont="BP", 
                           pvalueCutoff=0.05, qvalueCutoff=0.05)
)

# Compare UPREGULATED across all 3 contrasts
comp_up <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$up,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$up,
    Blum_vs_Water = sig_lists$blum_vs_water$up
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)

comp_down <- compareCluster(
  geneClusters = list(
    Water_vs_Ctrl = sig_lists$water_vs_ctrl$down,
    Blum_vs_Ctrl  = sig_lists$blum_vs_ctrl$down,
    Blum_vs_Water = sig_lists$blum_vs_water$down
  ),
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db, ont="BP", pvalueCutoff=0.05
)




dotplot(comp_up, showCategory=15, title = "Upregulated GO:BP")
dotplot(comp_down, showCategory=10, font.size=8, title = "Downregulated GO:BP")

# Individual contrast
dotplot(ego_up_all$water_vs_ctrl, showCategory=15, title = "Water vs Ctrl - Upregulated")
dotplot(ego_up_all$blum_vs_ctrl, showCategory=15, title = "Blum vs Ctrl - Upregulated")
dotplot(ego_up_all$blum_vs_water, showCategory=15, title = "Blum vs Water - Upregulated")

dotplot(ego_down_all$water_vs_ctrl, showCategory=15, title = "Water vs Ctrl - Downregulated")
dotplot(ego_down_all$blum_vs_ctrl, showCategory=15, title = "Blum vs Ctrl - Downregulated")
dotplot(ego_down_all$blum_vs_water, showCategory=15, title = "Blum vs Water - Downregulated")



