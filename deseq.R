library(tidyverse)
library(tximport)
library(DESeq2)
library(readr)

# Prepare your data
annot <- read_tsv("nv2.func-04.04.23.tsv")
tx2geneFile <- "Trinity.fasta.gene_trans_map"
tx2gene <- read.table(tx2geneFile, header = FALSE, col.names = c("gene_id", "transcript_id"))
tx2gene <- select(tx2gene, c(2,1))
samples <- list.dirs("kallisto", full.names = FALSE, recursive = FALSE)
files <- file.path("kallisto", samples, "abundance.tsv")
names(files) <- samples
sample_metadata <- data.frame(row.names = samples, condition = gsub("[0-9]+$", "", samples))

# import best blast hits and filter
bbh <- read_tsv("stomolophus_nv2.map.renamed.m8", col_names=FALSE)
names(bbh)<- c("query", "target", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen")
bbh <- bbh %>% filter(evalue <= 1e-2) %>% group_by(query) %>% arrange(desc(bits)) %>% 
       filter(row_number()==1) %>% ungroup()	

# Import transcript-level estimates and summarize to gene-level
# This preserves the fractional nature of the counts
txi <- tximport(files, 
                type = "kallisto", 
                tx2gene = tx2gene,  # transcript-to-gene mapping file
                countsFromAbundance = "no")  # Keeps original estimated counts

# Create DESeq2 dataset directly from tximport object
# DESeq2 internally handles these fractional counts appropriately
dds <- DESeqDataSetFromTximport(txi,
                                 colData = sample_metadata,
                                 design = ~ condition)

# Run DESeq2 - it knows how to work with these counts
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "SD", "SK"))
res <- as_tibble(res, rownames="gene") 

# join with best blast hits
res <- left_join(res, bbh, by=c("gene"="query"))
res <- res %>% filter(padj<=0.05) %>% left_join(annot %>% 
       select(geneID,gene_short_name,cdsName), by=c("target"="geneID"))

write_tsv(res, "DEG_results_all.csv")
