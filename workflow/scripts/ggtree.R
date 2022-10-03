
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(stringr)
library(ggplot2)
library(phytools)

if (interactive()) {
    setClass("snake", slots = list(input = "list", output = "list"))
    snakemake <- new("snake", input  = list(
            tree = "analysis/phylogeny/MCP_NCLDV_epa/epa_result.newick",
            fasta = "analysis/phylogeny/MCP_NCLDV.fasta",
            outgroups    = "metadata/queries/MCP_NCLDV_outgroups.faa",
            synonyms = "metadata/organisms.txt",
            hmm = Sys.glob("hmm_algae/*.hmm")
    ), output = list(
        image = "test.svg",
        jtree = "output/MCP_NCLDV.jtree"
    ))
}

with(snakemake@input, {
    tree_file     <<- tree
    fasta_file    <<- fasta
    synonyms_file <<- synonyms
    outgroup_file <<- outgroups
})
with(snakemake@output, {
    out_image_file <<- image
    out_jtree_file <<- jtree
})
with(snakemake@params, {
    outgroup_rooting <<- outgroup_rooting
})

read.fasta.headers <- function(fnames) {
    file.info(fnames) %>%
        filter(size > 0) %>%
        rownames %>%
        lapply(treeio::read.fasta) %>%
        lapply(names) %>%
        unlist %>%
        data.frame(title = .)
}

synonyms <- read.table(synonyms_file, header = T, sep = "\t", fill = T, na.strings = "") %>%
    mutate(Collapse = ifelse(is.na(Collapse), Name, Collapse))

headers <- read.fasta.headers(fasta_file) %>%
    extract(title, into = c("label", "ID"), regex = "^([^ ]+) ([^ ]+)", remove = F) %>%
    left_join(synonyms, by = "ID")

no_name <- filter(headers, is.na(Name)) %>%
    pull(label) %>%
    paste(collapse = ", ")
if (no_name != "") {
    print(paste("No aliases found for: ", no_name))
    quit(status = 1)
}

tree <- read.tree(tree_file)
tree <- phangorn::midpoint(tree, node.labels = "support")
if (outgroup_rooting) {
    outgroup_df <- read.fasta.headers(outgroup_file)
    outgroups <- with(outgroup_df, sub(" .*", "", title))
    tree <- ape::root(tree, node = MRCA(tree, outgroups), edgelabel = T, resolve.root = T)
}
tree <- as_tibble(tree) %>%
    mutate(support = ifelse(node %in% parent & label != "", label, NA)) %>%
    separate(support, into = c("SH_aLRT", "UFboot"), sep = "/", convert = T) %>%
    left_join(headers, by = "label") %>%
    mutate(label.show = Name) %>%
    mutate(isInternal = node %in% parent) %>%
    `class<-`(c("tbl_tree", "tbl_df", "tbl", "data.frame"))
tree_data <- as.treedata(tree)
write.jtree(tree_data, file = out_jtree_file)

ntaxa <- filter(tree, ! node %in% parent) %>% nrow

colors <- list(
    Haptophyta = "orange",
    Chlorophyta = "green",
    Streptophyta = "darkgreen",
    MAG = "purple",
    Stramenopiles = "brown",
    Cryptophyta = "red",
    Amoebozoa = "gold4",
    Euglenozoa = "yellow",
    Choanoflagellata = "darkslateblue",
    Glaucophyta = "cyan",
    Animals = "blue",
    Dinoflagellata = "gray50",
    Rhizaria = "gray30"
)

scaleClades <- function(p, df) {
    with(df, Reduce(function(.p, .node) {
        offs <- offspring(.p$data, .node)
        scale <- 0.5 / (nrow(offs) - 1)
        scaleClade(.p, .node, scale)
    }, node, p))
}
collapseClades <- function(p, df) {
    with(df, Reduce(function(.p, .node) {
        fill <- unlist(colors[Host[node == .node]])
        .p$data[.p$data$node == .node, "label.show"] <- label.show[node == .node]
        collapse(.p, .node, "mixed", fill = fill)
    }, node, p))
}
#labelClades <- function(p) {
#    with(df, Reduce(function(.p, .node) {
#        .p + geom_cladelab(node = .node, label = label[node == .node], align = T, offset = .2, textcolor = 'blue')
#    }, node, p))
#}

multi_species <- allDescendants(tree_data@phylo) %>%
    lapply(function(x) filter(tree, node %in% x)) %>%
    bind_rows(.id = "ancestor") %>%
    group_by(ancestor) %>%
    filter(n_distinct(Collapse, na.rm = T) == 1, sum(!isInternal) > 1) %>% # , !any(Group == "Haptophyta")) %>%
    ungroup %>%
    mutate(ancestor = as.numeric(ancestor)) %>%
    filter(! ancestor %in% node) %>%
    filter(!is.na(Collapse)) %>%
    group_by(ancestor, Collapse) %>%
    summarize(num_tips = sum(!isInternal), Host = first(na.omit(Host))) %>%
    mutate(label.show = sprintf("%s (%d)", Collapse, num_tips)) %>%
    rename(node = ancestor)
p <- ggtree(tree_data) +
    geom_nodepoint(aes(x = branch, subset = !is.na(UFboot) & UFboot >= 90, size = UFboot)) +
    geom_tiplab(aes(label = label.show), size = 4, align = T, linesize = 0) +
    geom_text2(aes(subset = node %in% multi_species$node, x = max(x, na.rm = T), label = label.show), nudge_x = 0.01, size = 4, hjust = 0) +
    geom_tippoint(aes(color = Host), size = 3) +
    geom_treescale(width = 0.5) +
    scale_size_continuous(limits = c(90, 100), range = c(1, 3)) +
    scale_shape_manual(values = seq(0,15)) +
    scale_color_manual(values = colors)

p <- scaleClades(p, multi_species)
p <- collapseClades(p, multi_species)
# p <- facet_plot(p, mapping = aes(x = as.numeric(as.factor(query.name)), shape = DESC), data = genes, geom = geom_point, panel = 'Genes')

ggsave(out_image_file, p, height = ntaxa * 0.1, width = 7, limitsize = F)

