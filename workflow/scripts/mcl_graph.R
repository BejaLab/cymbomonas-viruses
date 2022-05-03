
library(dplyr)
library(tidyr)
library(scatterpie)
library(tibble)
library(treeio)
library(igraph)
library(ggraph)
library(ape)

if (interactive()) {
    mcl_file <- "analysis/hh_mcl/mcl.txt"
    virus_metadata_file   <- "metadata/viruses.tsv"
    segment_metadata_file <- "metadata/viral_segments.tsv"
    segment_ffindex_files <- Sys.glob("analysis/PLV_segments/*.ffindex")
    virus_ffindex_files   <- Sys.glob("analysis/PLVs/*.ffindex")
    segment_pfam_files    <- Sys.glob("analysis/PLV_segments_hh/*-hhsearch-pfam.tsv")
    virus_pfam_files      <- Sys.glob("analysis/PLVs_hh/*-hhsearch-pfam.tsv")
    genes_cluster_file <- "metadata/genes_cluster.tsv"
    genes_pfam_file <- "metadata/genes_pfam.tsv"
    probab_threshold <- 80
    core_genes_file <- "pie.pdf"
    jtree_file      <- "output/MCP_PLV.jtree"
    segments_file   <- "metadata/viral_segments.tsv"
    reduced_tree_file <- "output/MCP_PLV_reduced.svg"
} else {
    with(snakemake@input, {
        mcl_file              <<- mcl
        virus_metadata_file   <<- virus_metadata
        segment_metadata_file <<- segment_metadata
        segment_pfam_files    <<- segment_pfam
        virus_pfam_files      <<- virus_pfam
        segment_ffindex_files <<- segment_ffindex
        virus_ffindex_files   <<- virus_ffindex
        genes_cluster_file    <<- genes_cluster
        genes_pfam_file       <<- genes_pfam
        jtree_file            <<- jtree
        segments_file         <<- segments
    })
    with(snakemake@params, {
        probab_threshold  <<- probab
    })
    with(snakemake@output, {
        data_file       <<- data
        core_genes_file <<- core_genes
        bipart_file     <<- bipartite
        reduced_tree_file <<- reduced_tree
    })
}
segments <- basename(segment_ffindex_files) %>% sub(".ffindex", "", .)
viruses  <- basename(virus_ffindex_files)   %>% sub(".ffindex", "", .)

genes_pfam <- read.table(genes_pfam_file, header = T, sep = "\t", fill = T, na.strings = "") %>%
    mutate(Comb = 1:n()) %>%
    separate_rows(Pfam, sep = ",") %>%
    group_by(Comb) %>%
    mutate(Comb.Num = n())
genes_cluster <- read.table(genes_cluster_file, header = T, sep = "\t")

metadata <- bind_rows(
    read.table(virus_metadata_file, sep = "\t", header = T),
    read.table(segment_metadata_file, sep = "\t", header = T)
) %>% `rownames<-`(.$short)
clusters <- readLines(mcl_file) %>%
    data.frame(ID = ., Cluster = 1:length(.)) %>%
    separate_rows(ID, sep = "\t") %>%
    left_join(genes_cluster, by = "ID") %>%
    group_by(Cluster) %>%
    mutate(Gene_Cluster = first(na.omit(Gene_Cluster)))

segment_pfam <- lapply(segment_pfam_files, read.table, header = T, sep = "\t", quote = "")
virus_pfam   <- lapply(virus_pfam_files, read.table, header = T, sep = "\t", quote = "")
pfam_all <- c(segment_pfam, virus_pfam) %>%
    bind_rows %>%
    select(-c(Q.query, T.ss_pred, T.hit, T.consensus, Q.consensus, T.ss_dssp, Q.ss_pred)) %>%
    extract(Hit.Description, into = "Hit.Name", regex = "; (.+?) ;", remove = F)
pfam <- filter(pfam_all, Probab >= probab_threshold)

to_gff <- function(.data) {
    mutate(.data, Source = "hhsearch", Feature = "domain", Strand = ".", Fname = ".", Attributes = paste(Hit.ID, Hit.Description)) %>%
        select(ID, Source, Feature, Q.Start, Q.End, Probab, Strand, Fname, Attributes)
}
most_freq <- function(.data) {
    na.omit(.data) %>%
        table(useNA = "always") %>%
        which.max %>%
        names
}
first_nona <- function(.data) {
    first(na.omit(.data))
}

segment_genes <- lapply(segment_ffindex_files, read.table, sep = "\t") %>%
    setNames(segments)
virus_genes   <- lapply(virus_ffindex_files, read.table, sep = "\t") %>%
    setNames(viruses)
data <- c(segment_genes, virus_genes) %>%
    bind_rows(.id = "Genome") %>%
    select(Genome, ID = 2) %>%
    left_join(pfam, by = c(ID = "Query")) %>%
    left_join(metadata, by = c(Genome = "short")) %>%
    left_join(genes_pfam, by = c(Hit.Name = "Pfam")) %>%
    group_by(Genome, ID, Comb) %>%
    mutate(Gene_Pfam = ifelse(n_distinct(Hit.ID) >= Comb.Num, Gene_Pfam, NA), Gene_Pfam = ifelse(is.na(No_threshold) | No <= No_threshold, Gene_Pfam, NA)) %>%
    left_join(clusters, by = "ID") %>%
    mutate(Cluster = as.character(ifelse(is.na(Cluster), ID, Cluster))) %>%
    group_by(Genome, ID, Cluster, clade, subclade, subsubclade) %>%
    summarize(Gene_Pfam = first_nona(Gene_Pfam), Gene_Cluster = first_nona(Gene_Cluster), Hit.ID = first(Hit.ID), Hit.Description = first(Hit.Description), Probab = first(Probab)) %>%
    group_by(Cluster) %>%
    mutate(Gene_Pfam = most_freq(Gene_Pfam), Gene_Cluster = most_freq(Gene_Cluster)) %>%
    mutate(Gene = ifelse(is.na(Gene_Cluster), Gene_Pfam, Gene_Cluster), Gene = ifelse(is.na(Gene), paste0("Cluster_", Cluster), Gene), Gene = ifelse(is.na(Gene), paste0("ID_", ID), Gene)) %>%
    ungroup
write.table(data, data_file, sep = "\t", quote = F, row.names = F)

virus_segments <- read.table(segments_file, sep = "\t", header = T)
tree <- read.jtree(jtree_file)
tip.data <- filter(tree@data, !isInternal) %>%
    mutate(label = sub(" .*", "", title)) %>%
    extract(title, into = c("scaffold", "gene_num", "genome", "start", "end"), regex = "^(.+?)_(\\d+) (.+?) \\[(\\d+) - (\\d+)\\]", remove = F, convert = T) %>%
    left_join(virus_segments, by = c("scaffold")) %>%
    filter(is.na(scaffold_start) | is.na(start) | start >= scaffold_start & end <= scaffold_end) %>%
    mutate(Name = ifelse(is.na(short), Name, short)) %>%
    {setNames(.$Name, .$label)}
tree@phylo$tip.label <- tip.data[tree@phylo$tip.label]
tips.to.keep <- filter(data, clade == "PgVV") %>%
    pull(Genome) %>%
    `[`(. %in% tree@phylo$tip.label) %>%
    c("PgVV")
tips.to.drop <- tree@phylo$tip.label[! tree@phylo$tip.label %in% tips.to.keep]
my.tree <- tree@phylo %>%
    keep.tip(tips.to.keep) %>%
    drop.tip(tips.to.drop) %>%
    ladderize %>%
    chronos(lambda = 0.1)
svg(reduced_tree_file)
plot(my.tree)
dev.off()
order.tips <- data.frame(tip.num = my.tree$edge[,2]) %>%
    filter(tip.num <= length(my.tree$tip.label)) %>%
    mutate(Genome = my.tree$tip.label[tip.num], num = n():1)

chosen_clades <- c("PgVV", "Endemic", "Mesomimi")
clade_levels <- c("Mesomimi", "Endemic", "Mavirus", "Virophage", "TVS", "PgVV")
ref_genomes <- c("Sputnik", "Mavirus_Spezl", "TVV_N1")

clades.data <- filter(data, clade %in% chosen_clades | Genome %in% ref_genomes) %>%
    mutate(Count = 1) %>%
    complete(Gene, nesting(Genome, clade, subclade, subsubclade), fill = list(Count = 0)) %>%
    group_by(Gene) %>%
    mutate(N_Subsubclades = n_distinct(subsubclade[clade == "PgVV" & Count > 0]), N_Genomes = n_distinct(Genome[clade == "PgVV" & Count > 0])) %>%
    filter(N_Subsubclades > 2) %>%
    mutate(clade = factor(clade, levels = clade_levels)) %>%
    left_join(order.tips, by = "Genome") %>%
    replace_na(list(num = 0)) %>%
    arrange(-N_Genomes, num, Gene, clade, desc(subclade), subsubclade, Genome) %>%
    ungroup %>%
    mutate(Gene = factor(Gene, levels = unique(Gene)), Genome = factor(Genome, levels = unique(Genome))) %>%
    filter(Count > 0)
ggplot_scatterpie <- function(.data, x.col, y.col, z.col, val.col, group.col) {
    x_uniq <- levels(.data[[x.col]])
    y_uniq <- levels(.data[[y.col]])
    .data <- ungroup(.data) %>%
        mutate(x = as.numeric(.[[x.col]]), y = as.numeric(.[[y.col]])) %>%
        rename(value = !!val.col)
    ggplot() +
        geom_scatterpie(aes(x = x, y = y), data = .data, long_format = T, cols = z.col, color = NA) +
        scale_x_continuous(breaks = 1:length(x_uniq), labels = x_uniq) +
        scale_y_continuous(breaks = 1:length(y_uniq), labels = y_uniq) +
        xlab(x.col) + ylab(y.col) #+
        #facet_grid(rows = group.col, scales = "free_y", space = "free_y")
}
p <- ggplot_scatterpie(clades.data, "Gene", "Genome", "Cluster", "Count", "subsubclade") +
    coord_equal() +
    theme_void() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text())
ggsave(core_genes_file, p, height = 12)

vertex_metadata <- list(
        Cluster = distinct(clades.data, Cluster, Gene) %>%
            mutate(label = ifelse(grepl("Cluster", Gene), NA, as.character(Gene)), subclade = "Cluster", subsubclade = "Cluster") %>%
            column_to_rownames("Cluster"),
        Genome = mutate(metadata, label = short)
    ) %>% bind_rows(.id = "Type")

g <- filter(data, clade == "PgVV") %>%
    distinct(Genome, Cluster) %>%
    group_by(Cluster) %>%
    filter(n() > 2) %>%
    mutate(present = T) %>%
    spread(Genome, present, F) %>%
    column_to_rownames("Cluster") %>%
    as.matrix %>%
    graph.incidence %>%
    set_vertex_attr("subclade",    value = vertex_metadata[V(.)$name,"subclade"]) %>%
    set_vertex_attr("subsubclade", value = vertex_metadata[V(.)$name,"subsubclade"]) %>%
    set_vertex_attr("label",       value = vertex_metadata[V(.)$name,"label"]) %>%
    set_vertex_attr("gene",        value = vertex_metadata[V(.)$name,"Gene"])

subsubclades <- list(
    blue = "#1F78B4",
    purple = "#6A3D9A",
    Isogal = "#B2DF8A",
    red = "#E31A1C",
    yellow = "yellow4",
    PleuPLV = "#33A02C",
    MLC = "#FDBF6F",
    SAF1 = "#FF7F00",
    RED2 = "#CAB2D6",
    Han1023 = "#A6CEE3",
    CCE = "#FB9A99",
    Delaware572 = "#B15928",
    Cluster = "black",
    "darkgray"
)

p <- ggraph(g, "fr") +
    geom_edge_link(alpha = .1) +
    geom_node_point(aes(shape = subclade, colour = subsubclade, size = type)) +
    geom_node_text(aes(filter = !type, label = label), size = 1.5, color = "red",   repel = T, nudge_y = -0.01, hjust = "right") +
    geom_node_text(aes(filter = type,  label = label), size = 2,   color = "black", repel = T, nudge_y = -0.01, hjust = "right") +
    scale_size_manual(values = c(1,2)) +
    scale_color_manual(values = subsubclades) +
    theme_graph(base_family = "sans")
ggsave(bipart_file, p, width = 10, height = 8)
