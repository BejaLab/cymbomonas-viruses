
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
            network  = "bak2/analysis/vcontact2/results/c1.ntw",
            genomes  = "bak2/analysis/vcontact2/results/genome_by_genome_overview.csv",
            profiles = "bak2/analysis/vcontact2/results/vConTACT_profiles.csv",
            segments = "metadata/viral_segments.tsv",
            viruses  = "metadata/viruses.tsv",
            colors   = "metadata/subclade_colors.txt"
        ),
        output = list(
            clustering = "bak2/output/vcontact2_clustering.svg"
        )
    )
}
with(snakemake@input, {
    network_file  <<- network
    genomes_file  <<- genomes
    profiles_file <<- profiles
    segments_file <<- segments
    viruses_file  <<- viruses
    colors_file   <<- colors
})

with(snakemake@output, {
    clustering_file <<- clustering
})

segments <- read.table(segments_file, sep = "\t", header = T)
viruses  <- bind_rows(
    read.table(segments_file, sep = "\t", header = T),
    read.table(viruses_file, sep = "\t", header = T)
)
subclades <- with(viruses, setNames(subclade, short))

sort_components <- function(.graph) {
    decompose(.graph) %>%
        `[`(order(sapply(., length), decreasing = T))
}

genomes <- read.csv(genomes_file, na.strings = "") %>%
    group_by(VC) %>%
    mutate(Cluster = ifelse(is.na(VC), Genome, paste(VC, paste(Genome, collapse = ", ")))) %>%
    left_join(viruses, by = c(Genome = "short"))
subclusters <- with(genomes, setNames(VC.Subcluster, Genome))

self_edges <- with(genomes, data.frame(A = Genome, B = Genome, weight = 1e-10))
network <- read.table(network_file, col.names = c("A","B","weight")) %>%
    bind_rows(self_edges) %>%
    mutate(weight = weight / 1000) %>%
    graph_from_data_frame(directed = F) %>%
    set_vertex_attr("subclade", value = subclades[V(.)$name]) %>%
    set_vertex_attr("subcluster", value = subclusters[V(.)$name])

colors <- read.table(colors_file, comment.char = "") %>%
    with(setNames(V2,V1)) %>%
    c("darkgray")

p <- ggraph(network, "fr") +
    geom_edge_link(aes(color = -weight, width = weight), alpha = .1) +
    geom_node_point(aes(color = subclade), size = 2) +
    geom_node_text(aes(label = ifelse(is.na(subcluster), name, sprintf("%s\n(%s)", name, subcluster))), repel = T, nudge_y = -0.01) +
    scale_color_manual(values = colors) # +
    #theme_graph()
ggsave(clustering_file, p, width = 10, height = 10)
