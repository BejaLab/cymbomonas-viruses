
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
            elements = "metadata/viral_elements.tsv",
            viruses  = "metadata/viruses.tsv"
        ),
        output = list(
            clustering = "output/vcontact2_clusrtring.svg"
        )
    )
}
with(snakemake@input, {
    network_file  <<- network
    genomes_file  <<- genomes
    profiles_file <<- profiles
    elements_file <<- elements
    viruses_file  <<- viruses
})

with(snakemake@output, {
    clustering_file <<- clustering
})

elements <- read.table(elements_file, sep = "\t", header = T)
viruses  <- bind_rows(
    read.table(elements_file, sep = "\t", header = T),
    read.table(viruses_file, sep = "\t", header = T)
)
clades <- with(viruses, setNames(clade, short))

sort_components <- function(.graph) {
    decompose(.graph) %>%
        `[`(order(sapply(., length), decreasing = T))
}

network <- read.table(network_file, col.names = c("A","B","weight")) %>%
    graph_from_data_frame(directed = F) %>%
    set_vertex_attr("clade", value = clades[V(.)$name])
p <- ggraph(network, "fr") +
    geom_edge_link(aes(color = -log2(weight), width = weight), alpha = .1) +
    geom_node_point(aes(colour = clade)) +
    geom_node_text(aes(label = name), repel = T, nudge_y = -0.01, hjust = "right") # +
    #theme_graph()
ggsave(clustering_file, p, width = 10, height = 10)

genomes <- read.csv(genomes_file, na.strings = "") %>%
    group_by(VC) %>%
    mutate(Cluster = ifelse(is.na(VC), Genome, paste(VC, paste(Genome, collapse = ", ")))) %>%
    left_join(viruses, by = c(Genome = "short")) %>%
    filter(! clade %in% c("Mesomimi", "Endemic"))

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
    set_vertex_attr("clade", value = clades[V(.)$name])

p <- ggraph(g, "stress") +
    geom_edge_link(alpha = .1) +
    geom_node_point(aes(colour = clade)) +
    geom_node_text(aes(filter = type, label = name), repel = T, nudge_y = -0.01, hjust = "right") +
    theme_graph()
# ggsave("plot2.pdf", p, width = 20, height = 20)

