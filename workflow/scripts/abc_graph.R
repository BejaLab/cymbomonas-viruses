
library(dplyr)
library(tidyr)

with(snakemake@input, {
    clu_tsv_file      <<- clu_tsv
    segment_tsv_files <<- segment_tsv
    virus_tsv_files   <<- virus_tsv
})
with(snakemake@params, {
    coverage_threshold <<- coverage
    probab_threshold   <<- probab
})
output_file <- unlist(snakemake@output)

clusters <- read.table(clu_tsv_file, sep = "\t", col.names = c("Cluster", "ID"))

segment_tsv <- lapply(segment_tsv_files, read.table, header = T, sep = "\t", quote = "")
virus_tsv   <- lapply(virus_tsv_files, read.table, header = T, sep = "\t", quote = "")
data <- c(segment_tsv, virus_tsv) %>%
    bind_rows %>%
    mutate(Q.Coverage = (Q.End - Q.Start + 1) / Q.Length * 100, T.Coverage = (T.End - T.Start + 1) / T.Length * 100) %>%
    filter(Q.Coverage >= coverage_threshold, T.Coverage >= coverage_threshold, Probab >= probab_threshold) %>%
    left_join(clusters, by = c(Hit.ID = "Cluster")) %>%
    select(Query, ID, Probab)
write.table(data, output_file, sep = "\t", row.names = F, col.names = F, quote = F)
