
library(dplyr)
library(tidyr)
library(ggplot2)

with(snakemake@input, {
    MCP_file <<- MCPs
    gff3_files <<- gff3
})
with(snakemake@params, {
    clades <<- clades
    complete_score <<- complete_score
    complete_lens <<- complete_lens
})
output_file <- unlist(snakemake@output)

lens_df <- data.frame(clade = names(complete_lens), len_threshold = unlist(complete_lens))

bed_cols <- c("scaffold", "offset", "end", "pident", "clade", "strand")
gff_cols <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

MCP <- read.table(MCP_file, col.names = bed_cols) %>%
    group_by(clade) %>%
    summarize(num_MCP = n())
gff3 <- lapply(gff3_files, read.table, col.names = gff_cols, sep = "\t") %>%
    setNames(clades) %>%
    bind_rows(.id = "clade") %>%
    mutate(len = end - start + 1) %>%
    left_join(lens_df, by = "clade") %>%
    mutate(is_complete = score > complete_score & ifelse(is.na(len_threshold), T, len >= len_threshold)) %>%
    group_by(clade) %>%
    summarize(total_len = sum(len), num_compl = sum(is_complete)) %>%
    left_join(MCP, by = "clade") %>%
    gather(key, value, -clade) %>%
    group_by(key) %>%
    mutate(pct = value / sum(value)) %>%
    mutate(ypos = 1 - cumsum(pct) + 0.5 * pct)

p <- ggplot(gff3, aes(x = key, y = pct, fill = clade)) +
    geom_col() +
    geom_text(aes(y = ypos, label = value), size=3, color = "white") +
    coord_polar(theta = "y") +
    theme_void()
ggsave(output_file, p)
