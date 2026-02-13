library(clusterProfiler)
library(org.Rn.eg.db)
library(DESeq2)

cm_results <- readRDS("DEG_Cardiomyocytes.rds")  # DESeqResults

# Pick one comparison; repeat later for others
res_water_ctrl <- cm_results$water_vs_ctrl   # class: DESeqResults (DataFrame)
res_blum_ctrl  <- cm_results$blum_vs_ctrl
res_blum_water <- cm_results$blum_vs_water

# Turn the chosen one into a data.frame
cm_df <- as.data.frame(res_water_ctrl)
cm_df$gene <- rownames(res_water_ctrl)  # gene symbols like Raet1e, Ulbp1, ...

id_map <- bitr(
  cm_df$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Rn.eg.db
)

cm_anno <- merge(cm_df, id_map, by.x = "gene", by.y = "SYMBOL")

sig_up   <- cm_anno$ENTREZID[cm_anno$padj < 0.05 & cm_anno$log2FoldChange >  1]
sig_down <- cm_anno$ENTREZID[cm_anno$padj < 0.05 & cm_anno$log2FoldChange < -1]
bg_genes <- unique(cm_anno$ENTREZID)

ego_up <- enrichGO(
  gene          = sig_up,
  universe      = bg_genes,
  OrgDb         = org.Rn.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(ego_up)  # Top BP terms upregulated
nrow(ego_up)  # How many significant?
dotplot(ego_up, showCategory=15)

ego_down <- enrichGO(
  gene          = sig_down,
  universe      = bg_genes,
  OrgDb         = org.Rn.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

head(ego_down)  # Top BP terms upregulated
nrow(ego_down)  # How many significant?
dotplot(ego_down, showCategory=15)
