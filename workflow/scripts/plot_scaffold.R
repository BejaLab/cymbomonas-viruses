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
})

clades <- unlist(snakemake@params["clades"])

output_file <- unlist(snakemake@output)

gff_cols <- c("seqnames", "source", "gff_type", "start", "end", "score", "strand", "frame", "comment")
tbl_cols <- c("gene", "target", "ID", "query", "full.sequence.E.value", "full.sequence.score", "full.sequence.bias", "best.1.domain.E.value", "best.1.domain.score", "best.1.domain.bias", "domain.number.estimation.exp", "domain.number.estimation.reg", "domain.number.estimation.clu", "domain.number.estimation.ov", "domain.number.estimation.env", "domain.number.estimation.dom", "domain.number.estimation.rep", "domain.number.estimation.inc", "description.of.target")

markers <- read_table(markers_file, col_names = c("gene", "marker", "threshold"))
gvogs <- read_table(gvog_file, comment = "#", col_names = tbl_cols) %>%
    left_join(markers, by = "gene")
gvog_markers <- filter(gvogs, full.sequence.score >= threshold) %>%
    select(ID, marker) %>%
    distinct(marker, .keep_all = T)

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
    mutate(track = recode(strand, `+` = 1, `-` = 2)) %>%
    replace_na(list(color = "#cccccc")) %>%
    mutate(shape = "rectangle", colour = color) %>%
    mutate(type = paste(COG, description, sep = ": ")) %>%
    mutate(label = "") %>%
    as_granges
marker_labels <- genes[!is.na(genes$marker)] %>%
    select(marker) %>%
    `$<-`("label", .$marker)

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
    feature_data = as_granges(genes), 
    feature_legend_title = "COG categories"
)

gmovizPlot(file_name = output_file, file_type = "svg", legends = legend, plotting_functions = {
    gmovizInitialise(
        genome,
        label_data = marker_labels,
        space_between_sectors = 25,
        start_degree = 90,
        sector_label_size = 1,
        track_height = 0.15,
        xaxis_spacing = 20000,
        xaxis_spacing_unit = "bp"
    )
    circos.par(track.margin = c(0, 0), cell.padding = c(-0.02, 0, -0.02, 0))
    drawFeatureTrack(c(genes, viruses), track_height = 0.1)
    drawLinegraphTrack(coverage_data, yaxis_increment = 1, ylim = c(0, max(coverage_data$coverage)+1), line_shade_colour = "#cccccc64", gridline_colour = "black")
})
