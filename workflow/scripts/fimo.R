library(dplyr)

input_file <- unlist(snakemake@input)
output_file <- unlist(snakemake@output)
with(snakemake@params, {
    motif_re <<- motif_re
    q_thresh <<- q_thresh
})
gff <- read.table(input_file, header = T) %>%
    filter(grepl(motif_re, matched_sequence), q.value < q_thresh) %>%
    mutate(description = sprintf("Name=%s;Alias=%s;pvalue=%f;qvalue=%f;sequence=%s;", motif_id, motif_alt_id, p.value, q.value, matched_sequence)) %>%
    mutate(source = "fimo", type = "nucleotide_motif", frame = ".") %>%
    select(sequence_name, source, type, start, stop, score, strand, frame, description)
write.table(gff, file = output_file, sep = "\t", row.names = F, col.names = F, quote = F)
