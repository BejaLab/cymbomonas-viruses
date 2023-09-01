library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

with(snakemake@input, {
    complement_gc_file <<- complement_gc
    clade_gc_files <<- clade_gc
    coverage_file <<- coverage
    scaffold_file <<- scaffold
})
with(snakemake@params, {
    clades <<- clades
    scaf_min <<- scaf_min
})
output_file <- unlist(snakemake@output)

clade_df <- lapply(clade_gc_files, read.table, col.names = c("scaffold", "size", "GC")) %>%
    setNames(clades) %>%
    bind_rows(.id = "clade") %>%
    group_by(scaffold) %>%
    mutate(longest = clade[which.max(size)]) %>%
    pivot_wider(id_cols = c(scaffold, longest), names_from = clade, values_from = c(size, GC))

complement_gc <- read.table(complement_gc_file, col.names = c("scaffold", "compl_size", "compl_gc"))

coverage <- read.table(coverage_file, col.names = c("scaffold", "scaf_size", "scaf_cov", "r"), fill = T) %>%
    filter(r > 0) %>%
    left_join(complement_gc, by = "scaffold") %>%
    replace_na(list(compl_size = 0)) %>%
    mutate(viral_size = scaf_size - compl_size) %>%
    mutate(viral_fract = viral_size / scaf_size) %>%
    left_join(clade_df, by = "scaffold") %>%
    group_by(longest) %>%
    mutate(med_scaf_cov = median(scaf_cov)) %>%
    ungroup %>%
    mutate(med_scaf_cov = first(med_scaf_cov[viral_fract == 0])) %>%
    mutate(scaf_cov_norm = scaf_cov / med_scaf_cov)

colors <- list(
    NCLDV = "#ccffaaff",
    PLV2  = "#5cc2eaff",
    PLVA  = "#d4865eff",
    PLVB  = "#fbbc56ff"
)

p <- ggplot() +
    geom_point(data = filter(coverage, viral_fract > 0, scaf_size > scaf_min), aes(x = viral_fract, y = scaf_cov_norm, size = scaf_size, color = longest)) +
    geom_violin(data = filter(coverage, viral_fract == 0), aes(x = 0, y = scaf_cov_norm), fill = NA) +
    xlim(c(0,1)) + ylim(c(0,5)) +
    scale_color_manual(values = colors) +
    xlab("Viral fraction in scaffold") + ylab("Normalized scaffold coverage") +
    theme_bw()
ggsave(output_file, p, height = 5, width = 7)
