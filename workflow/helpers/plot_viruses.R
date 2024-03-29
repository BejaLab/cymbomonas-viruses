library(dplyr)
library(ggplot2)
library(gggenomes)

genome_names <- c("jcf7180000153361_22007-6512","jcf7180000139581_28090-12920","AFUW01000088.1_22959-44618","AFUW01000076.1_27241-47092","jcf7180000142507_63534-39856","jcf7180000139874_69062-64535_48127-29878","jcf7180000145894_47802-25306","jcf7180000184515_27397-2173","jcf7180000150446_6083-27605","jcf7180000142501_16998-39223", "jcf7180000174485")
genome_files <- paste0("databases/manual_viruses/", genome_names, ".gbk")
gene_file <- "output/mcl_genes_60.tsv"

genes <- read.table(gene_file, header = T, sep = "\t") %>%
    filter(!grepl("Cluster_", Gene))
annotations <- lapply(genome_files, read_gbk) %>%
    bind_rows %>%
    left_join(genes, by = c(seq_id = "Genome", locus_tag = "ID")) %>%
    mutate(name = Gene)
cds <- filter(annotations, type == "CDS")

get_tirs <- function(df, left = T) {
    filter(df, type == "repeat_region", start < 100) %>%
        mutate(left = left, start = unlist(ifelse(left, start, data.frame(introns)[2,] + 1)), end = unlist(ifelse(left, data.frame(introns)[1,] - 1, end))) %>%
        mutate(strand = ifelse(left, "+", "-")) %>%
        select(-introns, -left)
}

tirs <- bind_rows(get_tirs(annotations, T), get_tirs(annotations, F))
p <- gggenomes(genes = annotations, feats = tirs) +
    geom_seq() +
    geom_feat(size = 8) +
    geom_gene(aes(fill = name)) +
    geom_gene_tag(aes(label = name), size = 2, nudge_y = 0.1, check_overlap = T) +    
    geom_seq_label() +
    theme(legend.position = "none")
ggsave("manual_viruses.pdf", p, height = 9, width = 3)
