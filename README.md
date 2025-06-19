# Tiwi_RNA_Seq_DESeq
Differential expression analysis for the Tiwi pathogenic hits


**Step 1: Create a Sample file to expression data using DESeq2**

      Rscript createSampleFileforDESeq.r /work/NagarajLab_BG/16_Tiwi_RNAseq/X201SC25050245-Z01-F002/03processing/star_rsem \
                                          /work/NagarajLab_BG/16_Tiwi_RNAseq/01_DifferntialExp \
                                          all_gene
