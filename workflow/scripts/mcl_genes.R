
library(dplyr)
library(tidyr)
library(scatterpie)
library(tibble)
library(treeio)
library(igraph)
library(ggraph)
library(ape)

with(snakemake@input, {
    mcl_file              <<- mcl
    virus_metadata_file   <<- virus_metadata
    virus_pfam_files      <<- virus_pfam
    virus_ffindex_files   <<- virus_ffindex
    genes_cluster_file    <<- genes_cluster
    genes_pfam_file       <<- genes_pfam
})
with(snakemake@params, {
    probab_threshold  <<- probab
})
with(snakemake@output, {
    data_file       <<- data
    bipart_file     <<- bipartite
})
viruses  <- basename(virus_ffindex_files)   %>% sub(".ffindex", "", .)

genes_pfam <- read.table(genes_pfam_file, header = T, sep = "\t", fill = T, na.strings = "") %>%
    mutate(Comb = 1:n()) %>%
    separate_rows(Pfam, sep = ",") %>%
    group_by(Comb) %>%
    mutate(Comb.Num = n())
genes_cluster <- read.table(genes_cluster_file, header = T, sep = "\t")

metadata <- read.table(virus_metadata_file, sep = "\t", header = T) %>%
    `rownames<-`(.$short)
clusters <- readLines(mcl_file) %>%
    data.frame(ID = ., Cluster = 1:length(.)) %>%
    separate_rows(ID, sep = "\t") %>%
    left_join(genes_cluster, by = "ID") %>%
    group_by(Cluster) %>%
    mutate(Gene_Cluster = first(na.omit(Gene_Cluster)))

pfam_all <- lapply(virus_pfam_files, read.table, header = T, sep = "\t", quote = "") %>%
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

data <- lapply(virus_ffindex_files, read.table, sep = "\t") %>%
    setNames(viruses) %>%
    bind_rows(.id = "Genome") %>%
    select(Genome, ID = 2) %>%
    left_join(pfam, by = c(ID = "Query")) %>%
    left_join(metadata, by = c(Genome = "short")) %>%
    left_join(genes_pfam, by = c(Hit.Name = "Pfam")) %>%
    group_by(Genome, ID, Comb) %>%
    mutate(Gene_Pfam = ifelse(n_distinct(Hit.ID) >= Comb.Num, Gene_Pfam, NA), Gene_Pfam = ifelse(is.na(No_threshold) | No <= No_threshold, Gene_Pfam, NA)) %>%
    left_join(clusters, by = "ID") %>%
    mutate(Cluster = as.character(ifelse(is.na(Cluster), ID, Cluster))) %>%
    ungroup %>%
    arrange(No, -Comb.Num) %>%
    group_by(Genome, ID, Cluster, clade) %>%
    summarize(Gene_Pfam = first_nona(Gene_Pfam), Gene_Cluster = first_nona(Gene_Cluster), Hit.ID = first(Hit.ID), Hit.Description = first(Hit.Description), Probab = first(Probab)) %>%
    group_by(Cluster) %>%
    mutate(Gene_Pfam = most_freq(Gene_Pfam), Gene_Cluster = most_freq(Gene_Cluster)) %>%
    mutate(Gene = ifelse(is.na(Gene_Cluster), Gene_Pfam, Gene_Cluster), Gene = ifelse(is.na(Gene), paste0("Cluster_", Cluster), Gene), Gene = ifelse(is.na(Gene), paste0("ID_", ID), Gene)) %>%
    ungroup
write.table(data, data_file, sep = "\t", quote = F, row.names = F)

vertex_metadata <- list(
        Cluster = distinct(data, Cluster, Gene) %>%
            mutate(label = ifelse(grepl("Cluster", Gene), NA, as.character(Gene))) %>%
            column_to_rownames("Cluster"),
        Genome = mutate(metadata, label = short)
    ) %>% bind_rows(.id = "Type")

g <- distinct(data, Genome, Cluster) %>%
    group_by(Cluster) %>%
    filter(n() > 2) %>%
    mutate(present = T) %>%
    spread(Genome, present, F) %>%
    column_to_rownames("Cluster") %>%
    as.matrix %>%
    graph.incidence %>%
    set_vertex_attr("label", value = vertex_metadata[match(V(.)$name, rownames(vertex_metadata)),"label"]) %>%
    set_vertex_attr("gene",     value = vertex_metadata[match(V(.)$name, rownames(vertex_metadata)),"Gene"]) %>%
    set_vertex_attr("node_type", value = vertex_metadata[match(V(.)$name, rownames(vertex_metadata)),"Type"])

set.seed(1234)
p <- ggraph(g, "fr") +
    geom_edge_link(alpha = .1) +
    geom_node_point(aes(shape = node_type, size = type)) +
    geom_node_text(aes(filter = !type, label = label), size = 1.5, color = "red",   repel = T, nudge_y = -0.01, hjust = "right") +
    geom_node_text(aes(filter = type,  label = label), size = 2,   color = "black", repel = T, nudge_y = -0.01, hjust = "right") +
    scale_size_manual(values = c(1,2)) +
    theme_graph(base_family = "sans")
ggsave(bipart_file, p, width = 10, height = 8)
