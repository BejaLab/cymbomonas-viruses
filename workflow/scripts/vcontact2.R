
library(dplyr)
library(tidyr)
library(tools)
library(ggraph)
library(ggforce)
library(igraph)
library(tidyverse)

if (interactive()) {
    setClass("snakemake", slots = list(input = "list", output = "list", params = "list", wildcards = "list"))
    snakemake <- new("snakemake",
        input  = list(
            network  = "analysis/vcontact2/results/c1.ntw",
            genomes  = "analysis/vcontact2/results/genome_by_genome_overview.csv",
            profiles = "analysis/vcontact2/results/vConTACT_profiles.csv",
            segments = "metadata/viral_segments.tsv",
            viruses  = "metadata/viruses.tsv"
        ),
        output = list(
            clustering = "output/vcontact2_clustering.svg"
        )
    )
}
with(snakemake@input, {
    network_file  <<- network
    genomes_file  <<- genomes
    profiles_file <<- profiles
    segments_file <<- segments
    viruses_file  <<- viruses
})

with(snakemake@output, {
    clustering_file <<- clustering
})

segments <- read.table(segments_file, sep = "\t", header = T)
viruses  <- bind_rows(
    read.table(segments_file, sep = "\t", header = T),
    read.table(viruses_file, sep = "\t", header = T)
) %>%
    mutate(subsubclade = ifelse(is.na(subsubclade), subclade, subsubclade))
subclades <- with(viruses, setNames(subclade, short))
subsubclades <- with(viruses, setNames(subsubclade, short))

sort_components <- function(.graph) {
    decompose(.graph) %>%
        `[`(order(sapply(., length), decreasing = T))
}

genomes <- read.csv(genomes_file, na.strings = "") %>%
    group_by(VC) %>%
    mutate(Cluster = ifelse(is.na(VC), Genome, paste(VC, paste(Genome, collapse = ", ")))) %>%
    left_join(viruses, by = c(Genome = "short")) %>%
    filter(! clade %in% c("Mesomimi", "Endemic"))
self_edges <- with(genomes, data.frame(A = Genome, B = Genome, weight = 0))
network <- read.table(network_file, col.names = c("A","B","weight")) %>%
    bind_rows(self_edges) %>%
    graph_from_data_frame(directed = F) %>%
    set_vertex_attr("subclade", value = subclades[V(.)$name]) %>%
    set_vertex_attr("subsubclade", value = subsubclades[V(.)$name])

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
    Endemic = "Darkgray",
    "Gray"
)
p <- ggraph(network, "fr") +
    geom_edge_link(aes(color = -weight, width = weight), alpha = .1) +
    geom_node_point(aes(shape = subclade, colour = subsubclade), size = 2) +
    geom_node_text(aes(label = name), repel = T, nudge_y = -0.01) + 
    scale_color_manual(values = subsubclades) # +
    #theme_graph()
ggsave(clustering_file, p, width = 10, height = 10)
q()

g <- read.csv(profiles_file) %>%
    filter(pc_id != "") %>%
    left_join(genomes, by = c(contig_id = "Genome")) %>%
    filter(!is.na(Cluster)) %>%
    distinct(contig_id, pc_id) %>%
    group_by(pc_id) %>%
    filter(n() > 1) %>%
    mutate(present = T) %>%
    spread(contig_id, present, F) %>%
    column_to_rownames("pc_id") %>%
    as.matrix %>%
    graph.incidence %>%
    set_vertex_attr("subclade",    value = subclades[V(.)$name]) %>%
    set_vertex_attr("subsubclade", value = subsubclades[V(.)$name])

p <- ggraph(g, "fr") +
    geom_edge_link(alpha = .1) +
    geom_node_point(aes(shape = subclade, colour = subsubclade)) +
    geom_node_text(aes(filter = type, label = name), repel = T, nudge_y = -0.01, hjust = "right") +
    theme_graph()
ggsave("plot2.pdf", p, width = 20, height = 20)

