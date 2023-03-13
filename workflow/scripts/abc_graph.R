
library(dplyr)
library(tidyr)

with(snakemake@input, {
    clu_tsv_file    <<- clu_tsv
    virus_tsv_files <<- virus_tsv
})
with(snakemake@params, {
    coverage <<- coverage
    probab   <<- probab
    coverage_q_frag <<- coverage_q_frag
    coverage_t_frag <<- coverage_t_frag
    identities_frag <<- identities_frag
    probab_frag     <<- probab_frag
})
output_file <- unlist(snakemake@output)

clusters <- read.table(clu_tsv_file, sep = "\t", col.names = c("Cluster", "ID"))

data <- lapply(virus_tsv_files, read.table, header = T, sep = "\t", quote = "") %>%
    bind_rows %>%
    mutate(Q.Coverage = (Q.End - Q.Start + 1) / Q.Length * 100, T.Coverage = (T.End - T.Start + 1) / T.Length * 100) %>%
    # filter(Probab >= probab, Q.Coverage >= coverage, T.Coverage >= coverage) %>%
    filter(Probab >= probab, Q.Coverage >= coverage & T.Coverage >= coverage | Q.Coverage >= coverage_q_frag & T.Coverage >= coverage_t_frag & Identities >= identities_frag & Probab >= probab_frag) %>%
    left_join(clusters, by = c(Hit.ID = "Cluster")) %>%
    select(Query, ID, Probab)
write.table(data, output_file, sep = "\t", row.names = F, col.names = F, quote = F)
