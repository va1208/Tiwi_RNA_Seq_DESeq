#### RNA Seq Pipeline - DESEQ2 ####
args <- commandArgs(trailingOnly = TRUE)
print(args)
##### import required packages #####
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
library("genefilter")
library("org.Hs.eg.db")
library("biomaRt")
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(pheatmap))
library(gplots)

#### function for the annotation ####
add.anns <- function(df, mart, ...) {
  nm <- rownames(df)
  anns <- getBM(
    attributes = c("ensembl_gene_id", "chromosome_name",
                   "start_position", "end_position", "strand", "external_gene_name",
                   "go_id", "name_1006", "definition_1006", "description"),
    filters = "ensembl_gene_id",
    values = nm,
    mart = mart
  )
  # Match order of annotations to df
  anns <- anns[match(nm, anns[, "ensembl_gene_id"]), ]

  colnames(anns) <- c("Ensembl_Gene_ID", "Chromosome", "Start", "End",
                      "Strand", "Gene_Symbol", "GO_IDs", "GO_Term", "GO_Definition", "Gene_Description")

  df <- cbind(anns, df[, 2:ncol(df)])
  rownames(df) <- nm
  return(df)
}

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
### create output directory ###
output_dir <- "/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp/deseq2_gene_results"
dir.create(output_dir, showWarnings = FALSE)


#define working directory #####
wd <- "/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp"
setwd(wd)

#### import the sample file #####
sample_table <- read.table("/work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp/all_gene_sample_file.tsv", header = T)
sample_table$condition <- as.factor(sample_table$condition)
sample_table$replicate <- as.factor(sample_table$replicate)

if(args[1] == 1){
#sub <- subset(sample_table, gene == "CLCN2")
# Import RSEM output
#txi <- tximport(sub$files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi <- tximport(sample_table$files, type = "rsem", txIn = FALSE, txOut = FALSE)
# Force all lengths to 1 to bypass the DESeq2 check
txi$length[txi$length <= 0 | is.na(txi$length)] <- 1

# Filter: remove genes with any 0 in abundance or length
#zero_or_unexpressed <- apply(txi$length == 0 | txi$abundance == 0, 1, any)
#table(zero_or_unexpressed)

#txi$length    <- txi$length[!zero_or_unexpressed, ]
#txi$abundance <- txi$abundance[!zero_or_unexpressed, ]
#txi$counts    <- txi$counts[!zero_or_unexpressed, ]

# Set sample names correctly
colnames(txi$counts) <- sample_table$sample
colnames(txi$abundance) <- sample_table$sample
colnames(txi$length) <- sample_table$sample
rownames(sample_table) <- sample_table$sample
sample_table$condition <- factor(sample_table$condition)

# Check: All lengths positive now?
#stopifnot(all(txi$length > 0))

### running all together to get the sample PC's to see how the samples are grouping together ######
### create dds object ###
dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

#log transformation
rld <- rlogTransformation(dds)
#heatmap

(mycols <- brewer.pal(6, "Dark2")[1:length(unique(sample_table$condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png(paste0(output_dir, "/qc-heatmap-samples.png"), w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[sample_table$condition], RowSideColors=mycols[sample_table$condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
#colData(rld)
#PCA custom
png(paste0(output_dir, "/PCA_samples.png"), w=2000, h=2000, res=300)
pcadata<-DESeq2::plotPCA(rld, intgroup=c("condition", "gene"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
ggplot(pcadata, aes(PC1, PC2, color=condition, shape=gene)) +
  geom_point(size=3) + theme_classic() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
#DESeq2::plotPCA(rld, intgroup="condition")
#dev.new()
saveRDS(dds, paste0(output_dir, "/dds_allSamples_together.RDS"))
saveRDS(rld, paste0(output_dir, "/rld_allSamples_together.RDS"))
rm(dds)
rm(rld)
} else {
### Do for gene by gene ###
gc()
gene_list <- unique(sample_table$gene)

for (gene_name in gene_list) {
  message("Processing gene: ", gene_name)
  
  sub <- subset(sample_table, gene == gene_name)
  print(sub)
  sub$condition <- factor(sub$condition)
  #sub$condition <- relevel(factor(sub$condition), ref = "WT")
  rownames(sub) <- sub$sample

# Set up valid comparisons
  gene_conditions <- levels(droplevels(sub$condition))
  print(gene_conditions)
  valid_pairs <- list()

if (gene_name != "STXBP2" && all(c("WT", "MUT") %in% gene_conditions)) {
    valid_pairs <- list(c("MUT", "WT"))
  } else if (gene_name == "STXBP2") {
    if (all(c("WT", "MUT1") %in% gene_conditions)) valid_pairs <- append(valid_pairs, list(c("MUT1", "WT")))
    if (all(c("WT", "MUT2") %in% gene_conditions)) valid_pairs <- append(valid_pairs, list(c("MUT2", "WT")))
  }

if (length(valid_pairs) == 0) {
    message("Skipping gene ", gene_name, ": no valid condition pairs.")
    next
  }
txi <- tximport(sub$files, type = "rsem", txIn = FALSE, txOut = FALSE)
# Force all lengths to 1 to bypass the DESeq2 check
txi$length[txi$length <= 0 | is.na(txi$length)] <- 1
gc()
# Run DESeq2 once for all conditions
dds <- DESeqDataSetFromTximport(txi, colData = sub, design = ~condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
if(gene_name == "STXBP2") {
  dds$condition <- factor(dds$condition, levels = c("WT", "MUT1", "MUT2"))
} else {
  dds$condition <- factor(dds$condition, levels = c("WT", "MUT"))
}
dds$condition <- relevel(dds$condition, ref = "WT")  ## set WT as reference ####
### run DESeq2
dds <- DESeq(dds)

gc()
raw_counts <- counts(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  write.table(raw_counts,
              file = file.path(output_dir, paste0(gene_name, "_counts_raw.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(norm_counts,
              file = file.path(output_dir, paste0(gene_name, "_counts_normalized.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)

print(table(dds$condition))
print(resultsNames(dds))
# Loop through valid comparisons
  for (pair in valid_pairs) {
    mut <- pair[1]
    wt <- pair[2]
    contrast_label <- paste0(mut, "_vs_", wt)
    message("  Comparing: ", contrast_label)
    
    ## Get DE results ###
    res <- results(dds, contrast = c("condition", mut, wt))
    print(resultsNames(dds))
    # Shrink LFC ###
    coef_name <- paste0("condition_", mut, "_vs_", wt)
    print(coef_name)
    res_shrink <- lfcShrink(dds, coef = coef_name, type = "ashr")
    out <- add.anns(res_shrink, mart=ensembl)
    write.table(as.data.frame(res),
                file = file.path(output_dir, paste0(gene_name, "_", contrast_label, "_deseq2_results_raw.tsv")),
                sep = "\t", quote = FALSE)

    # Save shrunken results
    write.table(as.data.frame(out),
                file = file.path(output_dir, paste0(gene_name, "_", contrast_label, "_deseq2_results_shrunken.tsv")),
                sep = "\t", quote = FALSE)
  }
}
