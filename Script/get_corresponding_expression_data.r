### get expression data ###
library(data.table)
library(dplyr)

#### import the id table ####
id_table <- fread("/mnt/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp/id_gene.txt")

### raw counts ###
raw_counts_files <- list.files(pattern = "*_counts_raw.tsv", full.names = TRUE, path = "/mnt/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp/deseq2_gene_results/")
raw_count_names <- gsub("_counts_raw.tsv", "", basename(raw_counts_files))
raw_counts <- lapply(raw_counts_files, fread)
names(raw_counts) <- raw_count_names

### import normalized counts ###
norm_counts_files <- list.files(pattern = "*_counts_normalized.tsv", full.names = TRUE, path = "/mnt/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp/deseq2_gene_results/")
norm_counts_names <- gsub("_counts_normalized.tsv", "", basename(norm_counts_files))
norm_counts <- lapply(norm_counts_files, fread)
names(norm_counts) <- norm_counts_names

### subset the raw counts and normalized counts based on the id_table ###
raw_counts_filter <- list()
for (gene in names(raw_counts)) {
    sub <- subset(id_table, Gene == gene)
    out <- raw_counts[[gene]] %>% filter(V1 %in% sub$ensembl_id)
    out$V1 <- paste0(out$V1, "_raw")
    out$Gene <- gene
    raw_counts_filter[[gene]] <- out
}

norm_counts_filter <- list()
for (gene in names(norm_counts)) {
    sub <- subset(id_table, Gene == gene)
    out <- norm_counts[[gene]] %>% filter(V1 %in% sub$ensembl_id)
    out$V1 <- paste0(out$V1, "_norm")
    out$Gene <- gene
    norm_counts_filter[[gene]] <- out
}

#### combine the both for each gene ####
combined_counts <- list()
for (gene in names(raw_counts_filter)) {
    combined <- rbind(raw_counts_filter[[gene]], norm_counts_filter[[gene]])
    combined_counts[[gene]] <- combined
}

### save a xlsx file with all the counts ###
lapply(combined_counts, function(x) {
    write.table(x, file = paste0("/mnt/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp/deseq2_gene_results/", x$Gene[1], "_counts_raw_norm.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
})