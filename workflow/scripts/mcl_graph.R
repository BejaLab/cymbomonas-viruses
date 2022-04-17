
library(dplyr)
library(tidyr)

with(snakemake@input, {
    mcl_file              <<- mcl
    segment_pfam_files    <<- segment_pfam
    virus_pfam_files      <<- virus_pfam
    segment_ffindex_files <<- segment_ffindex
    virus_ffindex_files   <<- virus_ffindex
})
with(snakemake@params, {
    probab_threshold  <<- probab
})
with(snakemake@wildcards, {
    segments <<- segment
    viruses  <<- virus
})

output_file <- unlist(snakemake@output)

clusters <- readLines(mcl_file) %>%
    data.frame(ID = ., Cluster = 1:length(.)) %>%
    separate_rows(ID, sep = "\t")

segment_pfam <- lapply(segment_pfam_files, read.table, header = T, sep = "\t", quote = "")
virus_pfam   <- lapply(virus_pfam_files, read.table, header = T, sep = "\t", quote = "")
pfam <- c(segment_pfam, virus_pfam) %>%
    bind_rows %>%
    select(-c(Q.query, T.ss_pred, T.hit, T.consensus, Q.consensus, T.ss_dssp)) %>%
    filter(Probab >= probab_threshold)

segment_genes <- lapply(segment_ffindex_files, read.table, sep = "\t") %>%
    setNames(segments)
virus_genes   <- lapply(virus_ffindex_files, read.table, sep = "\t") %>%
    setNames(viruses)
data <- c(segment_genes, virus_genes) %>%
    bind_rows(.id = "Genome") %>%
    select(Genome, ID = 2) %>%
    left_join(clusters, by = "ID") %>%
    left_join(pfam, by = c(ID = "Query"))
write.table(data, output_file, sep = "\t", row.names = F, col.names = F, quote = F)
