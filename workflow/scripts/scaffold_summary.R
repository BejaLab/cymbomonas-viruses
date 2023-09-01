library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

with(snakemake@input, {
    multi_exon_file <- "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.multiexon.tab"
    complement_gc_file <- "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed.gc"
    clade_gc_files <- Sys.glob("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-*.gff3.gc")
    coverage_file <- "coverage/PLY262_illumina_001_v_PLY262_masurca_assembly.bam.coverage"
    MCP_file <- "analysis/locate_viruses/all_MCPs/cymbo_MaSURCA_assembly_2020.bed"
    gff3_files <- Sys.glob("analysis/locate_viruses/segments/*-*.gff3")
})
with(snakemake@params, {
    clades <<- clades
})
output_file <- unlist(snakemake@output)

bed_cols <- c("scaffold", "offset", "end", "pident", "clade", "strand")
MCP_scaf <- read.table(MCP_file, col.names = bed_cols) %>%
    mutate(MCP_loc = sprintf("%s:%d-%d(%s)", clade, offset + 1, end, strand)) %>%
    group_by(scaffold, clade) %>%
    summarize(MCP = n())

multi_exon <- read.table(multi_exon_file, sep = ",", col.names = c("scaffold", "start", "end", "comment")) %>%
    extract(comment, into = "product", regex = "product=([^;]+)", remove = F) %>%
    extract(comment, into = "dbxref", regex = "Dbxref=([^;]+)", remove = F) %>%
    extract(comment, into = "note", regex = "Note=([^;]+)", remove = F) %>%
    extract(comment, into = "ID", regex = "ID=([^;]+)") %>%
    replace_na(list(dbxref = "", note = "")) %>%
    mutate(dbxref = paste(dbxref, note, sep = ",")) %>%
    mutate(dbxref = gsub("^,|,$|", "", dbxref)) %>%
    filter(!grepl("hypothetical protein", product) | dbxref != "") %>%
    mutate(product = paste(ID, product, ifelse(dbxref == "", "", sprintf("(%s)", dbxref)))) %>%
    group_by(scaffold) %>%
    summarize(products = paste(product, collapse = "; "))

complement_gc <- read.table(complement_gc_file, col.names = c("scaffold", "compl_size", "compl_gc"))

coverage <- read.table(coverage_file, col.names = c("scaffold", "scaf_size", "scaf_cov", "r"), fill = T) %>%
    left_join(complement_gc, by = "scaffold") %>%
    replace_na(list(compl_size = 0)) %>%
    mutate(viral_size = scaf_size - compl_size) %>%
    mutate(viral_fract = viral_size / scaf_size)

clade_df <- lapply(clade_gc_files, read.table, col.names = c("scaffold", "size", "GC")) %>%
    setNames(clades) %>%
    bind_rows(.id = "clade") %>%
    left_join(MCP_scaf, by = c("scaffold", "clade")) %>%
    group_by(scaffold) %>%
    mutate(longest = clade[which.max(size)]) %>%
    pivot_wider(id_cols = c(scaffold, longest), names_from = clade, values_from = c(size, GC, MCP)) %>%
    left_join(coverage, by = "scaffold")

coverage_df <- left_join(coverage, select(clade_df, scaffold, longest), by = "scaffold") %>%
    group_by(longest) %>%
    mutate(avg_scaf_cov = mean(scaf_cov, trim = 0.1)) %>%
    ungroup %>%
    mutate(avg_scaf_cov = first(avg_scaf_cov[viral_fract == 0])) %>%
    mutate(scaf_cov_norm = scaf_cov / avg_scaf_cov) %>%
    filter(scaf_size > 10000)

colors <- list(
    NCLDV = "#ccffaaff",
    PLV2  = "#5cc2eaff",
    PLVA  = "#d4865eff",
    PLVB  = "#fbbc56ff"
)

p <- ggplot() +
    geom_point(data = filter(coverage_df, viral_fract > 0), aes(x = viral_fract, y = scaf_cov_norm, size = scaf_size, color = longest)) +
    geom_violin(data = filter(coverage_df, viral_fract == 0), aes(x = 0, y = scaf_cov_norm), fill = NA) +
    xlim(c(0,1)) + ylim(c(0,5)) +
    scale_color_manual(values = colors) +
    theme_bw()
ggsave("tmp.pdf", p)

p <- ggplot(clade_cov, aes(x = clade, y = scaf_cov, color = clade)) +
    geom_violin()
    # geom_jitter(width = 0.15)

per_pos <- read.table(pipe("grep jcf7180000139292 coverage/depth.PLY262_illumina_001_v_PLY262_masurca_assembly.bam.out"), col.names=c("scaffold","pos","cov"))

p <- ggplot(per_pos) +
    geom_line(aes(x = pos, y = rollmean(cov, 1000, na.pad = TRUE))) +
    ylab("coverage") + xlab("position")

ggsave("per_pos.png", p)
