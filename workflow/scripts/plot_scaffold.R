library(gmoviz)
library(plyranges)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

options(readr.show_col_types = F)

with(snakemake@input, {
    eggnog_file <<- eggnog
    genes_file <<- genes
    cog_file <<- cogs
    fna_file <<- fna
    cov_file <<- cov
    cov_total_file <<- cov_total
    virus_files <<- viruses
    colors_file <<- colors
    gvog_file <<- gvog
    markers_file <<- markers
    custom_domains_file <<- custom_domains
    iprscan_file <<- iprscan
})

clades <- unlist(snakemake@params["clades"])

output_file <- unlist(snakemake@output)

gff_cols <- c("seqnames", "source", "gff_type", "start", "end", "score", "strand", "frame", "comment")

read_hmmer <- function(fname) {
    col.names <- c("gene", "target", "ID", "query", "full.sequence.E.value", "full.sequence.score", "full.sequence.bias", "best.1.domain.E.value", "best.1.domain.score", "best.1.domain.bias", "domain.number.estimation.exp", "domain.number.estimation.reg", "domain.number.estimation.clu", "domain.number.estimation.ov", "domain.number.estimation.env", "domain.number.estimation.dom", "domain.number.estimation.rep", "domain.number.estimation.inc", "description.of.target")
    numeric.cols <- which(col.names == "full.sequence.E.value")        : which(col.names == "domain.number.estimation.exp")
    integer.cols <- which(col.names == "domain.number.estimation.reg") : which(col.names == "domain.number.estimation.inc")
    readLines(fname) %>%
        data.frame(line = .) %>%
        filter(!grepl("^#", line)) %>%
        separate(line, into = col.names, sep = " +", extra = "merge", convert = F) %>%
        mutate_at(numeric.cols, as.numeric) %>%
        mutate_at(integer.cols, as.integer)
}

custom_domains <- read_tsv(custom_domains_file, col_names = c("signature", "alias"))
ipr_scan <- read_tsv(iprscan_file, col_names = c("ID", "checksum", "prot_len", "db", "target", "definition", "start", "end", "score", "flag", "date", "ipr", "description")) %>%
    left_join(custom_domains, by = c(ipr = "signature")) %>%
    left_join(custom_domains, by = c(target = "signature")) %>%
    mutate(alias = ifelse(is.na(alias.x), alias.y, alias.x)) %>%
    filter(!is.na(alias)) %>%
    distinct(ID, alias) %>%
    group_by(ID) %>%
    summarize(alias = paste(alias, collapse = "+"))

markers <- read_table(markers_file, col_names = c("gene", "marker", "threshold"))
gvogs <- read_hmmer(gvog_file) %>%
    left_join(markers, by = "gene")

gvog_markers <- filter(gvogs, full.sequence.score >= threshold) %>%
    select(ID, marker)

total_cov <- read_tsv(cov_total_file, col_names = c("scaffold", "size", "coverage", "fraction")) %>%
    filter(fraction > 0)

cogs <- read_tsv(cog_file)
eggnog <- read_tsv(eggnog_file, comment = "##", na = "-")

genes <- read_tsv(genes_file, col_names = gff_cols, comment = "#") %>%
    filter(gff_type == "transcript") %>%
    select(-score, -frame) %>%
    extract(comment, into = "ID", regex = "ID=([^;]+)") %>%
    left_join(eggnog, by = c(ID = "#query")) %>%
    extract(COG_category, into = "COG", regex = "(.)", remove = F) %>%
    left_join(cogs, by = c(COG = "COG")) %>%
    left_join(gvog_markers, by = "ID") %>%
    left_join(ipr_scan, by = "ID") %>%
    mutate(track = recode(strand, `+` = 1, `-` = 2)) %>%
    replace_na(list(color = "#cccccc")) %>%
    mutate(colour = color) %>%
    mutate(shape = case_when(alias == "RT" & strand == "+" ~ "forward_arrow", alias == "RT" & strand == "-" ~ "reverse_arrow", T ~ "rectangle")) %>%
    mutate(type = paste(COG, description, sep = ": ")) %>%
    mutate(label = "") %>%
    as_granges

gene_labels <- genes[!is.na(genes$marker)] %>%
    select(marker) %>%
    `$<-`("label", .$marker)
if (length(gene_labels) > 0) {
    gene_labels$color <- "blue"
} else {
    gene_labels <- genes[!is.na(genes$alias)] %>%
        select(alias) %>%
        `$<-`("label", .$alias)
    gene_labels$color <- "red"
}

window_size <- 1000
coverage_data <- read.table(cov_file, col.names = c("seqnames","pos","coverage")) %>%
     mutate(bin = cut_width(pos, window_size)) %>%
     group_by(seqnames, bin) %>%
     summarize(coverage = mean(coverage) / median(total_cov$coverage), start = mean(pos), end = mean(pos), .groups = "drop") %>%
     select(-bin) %>%
     as_granges

genome <- getIdeogramData(fasta_file = fna_file)

clade_colors <- read_tsv(colors_file, col_names = c("label", "colour"))
viruses <- lapply(virus_files, read_tsv, col_names = gff_cols, comment = "#") %>%
    setNames(clades) %>%
    bind_rows(.id = "label") %>%
    filter(seqnames == genome@seqinfo@seqnames) %>%
    mutate(shape = recode(strand, `+` = "forward_arrow", `-` = "reverse_arrow")) %>%
    mutate(track = 3) %>%
    left_join(clade_colors, by = "label") %>%
    mutate(width = end - start + 1) %>%
    group_by(label) %>%
    filter(sum(width) / genome@ranges@width < 0.5) %>%
    as_granges
legend <- makeLegends(
    feature_legend = T,
    feature_data = arrange(genes, COG),
    feature_legend_title = "COG categories"
)

gmovizPlot(file_name = output_file, file_type = "svg", legends = legend, plotting_functions = {
    gmovizInitialise(
        genome,
        label_data = gene_labels,
        label_colour = gene_labels$color,
        space_between_sectors = 25,
        start_degree = 90,
        sector_label_size = 1,
        track_height = 0.15,
        xaxis_spacing = 20000,
        xaxis_spacing_unit = "bp"
    )
    circos.par(track.margin = c(0, 0), cell.padding = c(-0.02, 0, -0.02, 0))
    drawFeatureTrack(c(genes, viruses), track_height = 0.1)
    drawLinegraphTrack(coverage_data, yaxis_increment = 1, ylim = c(0, max(coverage_data$coverage) + 0.5), line_shade_colour = "#cccccc64", gridline_colour = "black")
})
