### import required packages ###
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

#### set working directory ###
wd <- "/mnt/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp/deseq2_gene_results"
setwd(wd)

### create a directory for pathway analysis ###
if (!dir.exists("pathway_analysis")) {
  dir.create("pathway_analysis")
}

#### import the data ###
files <- list.files(path = wd, pattern = "_deseq2_results_shrunken.tsv", full.names = TRUE)
data_list <- lapply(files, fread, sep = "\t", stringsAsFactors = FALSE)
names(data_list) <- gsub("_deseq2_results_shrunken.tsv", "", basename(files))   

#### subset the data for GO analysis ####
go_data_list <- lapply(data_list, function(df) {
  df <- df %>%
    filter(!is.na(log2FoldChange) & !is.na(padj)) %>%
    filter(padj < 0.05) %>%
    dplyr::select(Gene_Symbol, Ensembl_Gene_ID, log2FoldChange, padj)
  
  # Convert to data.table for compatibility with clusterProfiler
  setDT(df)
  return(df)
})

### check the dimensions of the data ###
lapply(go_data_list, dim)

### remove the data if it has less than 50 genes ###
go_data_list <- go_data_list[sapply(go_data_list, nrow) >= 50]
final_data_list <- names(go_data_list)

### perform Over representation analysis ###
go_list <- list()
for (i in seq_along(go_data_list)) {
  sub_fc <- go_data_list[[i]]
  gene_list <- as.vector(sub_fc$Ensembl_Gene_ID)
  message(paste("Running for", final_data_list[i]))
  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = gene_list,
                                      columns = "ENTREZID",
                                      keytype = "ENSEMBL")
  colnames(entrez_ids) <- c("Ensembl_Gene_ID", "ENTREZID")
  
  sub_fc <- merge(sub_fc, entrez_ids, by = "Ensembl_Gene_ID", all.x = TRUE)
  sub_fc <- na.omit(sub_fc)
  final_list <- setNames(sub_fc$log2FoldChange, sub_fc$ENTREZID)
  gene_final <- names(final_list)[abs(final_list) > 1]
  message(paste0("Number of genes for ", final_data_list[i], ": ", length(gene_final)))
  ggo <- enrichGO(gene = gene_final,
                  OrgDb = org.Hs.eg.db,
                  ont = "all",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  saveRDS(ggo, paste0("pathway_analysis/ORA_go_result_", final_data_list[i], ".RDS"))
  go_list[[final_data_list[i]]] <- ggo
  write.csv(ggo@result, paste0("pathway_analysis/ORA_go_result_", final_data_list[i], ".csv"), quote = FALSE, row.names = FALSE)
}

### create a plot for the each go analysis ###
for (i in seq_along(go_list)) {
    go <- go_list[[i]]
    p1 <- dotplot(go, showCategory = 30) + 
      ggtitle(paste("ORA Analysis for", final_data_list[i])) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  png(paste0("pathway_analysis/ORA_Dotplot_", final_data_list[i], ".png"), width = 2500, height = 5000, res = 300, units = "px")
  print(p1)
  dev.off() 
}

#### run gene enrichment analyss ### gseGO
gsea_list <- list() 
for (i in seq_along(go_data_list)) {
  sub_fc <- go_data_list[[i]]
  gene_list <- as.vector(sub_fc$Ensembl_Gene_ID)

  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = gene_list,
                                      columns = "ENTREZID",
                                      keytype = "ENSEMBL")
  colnames(entrez_ids) <- c("Ensembl_Gene_ID", "ENTREZID")
  
  sub_fc <- merge(sub_fc, entrez_ids, by = "Ensembl_Gene_ID", all.x = TRUE)
  sub_fc <- na.omit(sub_fc)
  
  gene_list_named <- setNames(sub_fc$log2FoldChange, sub_fc$ENTREZID)
  gene_list_named <- sort(gene_list_named, decreasing = TRUE)

  gsea_go <- gseGO(geneList = gene_list_named,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   keyType = "ENTREZID",
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = TRUE)
  
  saveRDS(gsea_go, paste0("pathway_analysis/GSEa_go_result_", final_data_list[i], ".RDS"))
  gsea_list[[final_data_list[i]]] <- gsea_go  
  write.csv(gsea_go@result, paste0("pathway_analysis/GSEA_GO_result_", final_data_list[i], ".csv"), quote = FALSE, row.names = FALSE)
}

### create a plot for the each go analysis ###
for (i in seq_along(gsea_list)) {
    go <- gsea_list[[i]]
    p1 <- dotplot(go, showCategory = 30) + 
      ggtitle(paste("GSEA Analysis for", final_data_list[i])) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  png(paste0("pathway_analysis/GSEA_Dotplot_", final_data_list[i], ".png"), width = 2500, height = 5000, res = 300, units = "px")
  print(p1)
  dev.off() 
}

