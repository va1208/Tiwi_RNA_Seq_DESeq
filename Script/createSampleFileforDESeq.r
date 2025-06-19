commandArgs(trailingOnly = TRUE)
print(args)

#### explanation for the args ####
### args 1  - file path for star_rsem ###
### args 2 - path to save the sample file ### dont put last slash in the path
### args 3 - output prefix
#### create a sample data frame for the files ####

files <- list.files(path = args[1], pattern="*.genes.results", full.names=TRUE) ### get the file path and file name ###
legnth(files)

### get sample name from the file names ###
sample_names <- gsub(".genes.results", "", basename(files))
length(sample_names)

print(length(files) == length(sample_names)) ## sanity check 

### create condition from the file names ###
condition <- ifelse(grepl("WT", sample_names), "WT",
              ifelse(grepl("MUT1", sample_names), "MUT1",
              ifelse(grepl("MUT2", sample_names), "MUT2", "MUT")))


print(length(files) == length(condition)) ## sanity check

### create dataframe 
samples <- data.frame(sample = sample_names, 
                            condtion = condition,
                            files = files,
                            stringsAsFactors = FALSE)

print(dim(samples))

write.table(samples, paste0(args[2], "/", args[3] "_sample_file.tsv"), quote = F, row.names = F)
print("all gene sample file saved in the given folder")

### spilit by gene ###
gene_names <- sub("_.*", "", sample_names)

samples <- data.frame(gene = gene_names,
                            sample = sample_names, 
                            condtion = condition,
                            files = files,
                            stringsAsFactors = FALSE)

### create a folder to save the file #####
if (!dir.exists("sample_metadata_by_gene")) {
  dir.create("sample_metadata_by_gene")
}

### get unique gene####
unique_genes <- unique(samples_df$gene)

for (g in unique_genes) {
  gene_df <- samples_df[samples_df$gene == g, ]
  output_file <- paste0(args[2], "/","sample_metadata_by_gene", "/", g, "_samples.tsv"))
  write.table(gene_df, output_file, row.names = F, quote = F)
  message("Written metadata for gene: ", g)
}
