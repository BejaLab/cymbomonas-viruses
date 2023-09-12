library(dplyr)
library(tidyr)

with(snakemake@input, {
    mRNA_file <<- mRNA # "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.mRNA.tab"
    complement_gc_file <<- complement_gc # "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed.gc"
    clade_gc_files <<- clade_gc # Sys.glob("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-*.gff3.gc")
    coverage_file <<- coverage # "coverage/PLY262_illumina_001_v_PLY262_masurca_assembly.bam.coverage"
    MCP_file <<- MCP # "analysis/locate_viruses/all_MCPs/cymbo_MaSURCA_assembly_2020.bed"
    gff3_files <<- gff3 # Sys.glob("analysis/locate_viruses/segments/*-*.gff3")
    proteinortho_file <<- proteinortho # "../genome_annotation/analysis/proteinortho/proteinortho.proteinortho.tsv"
})
with(snakemake@params, {
    clades <<- clades
    min_species <<- min_species
})
output_file <- unlist(snakemake@output)

bed_cols <- c("scaffold", "offset", "end", "pident", "clade", "strand")
MCP_scaf <- read.table(MCP_file, col.names = bed_cols) %>%
    mutate(MCP_loc = sprintf("%s:%d-%d(%s)", clade, offset + 1, end, strand)) %>%
    group_by(scaffold, clade) %>%
    summarize(MCP = n())

proteinortho <- read.table(proteinortho_file, sep = "\t", comment.char = "", na.strings = "*", header = T) %>%
    select(protein_id = annot.faa, num_species = X..Species, num_genes = Genes, alg_conn = Alg..Conn.) %>%
    filter(!is.na(protein_id)) %>%
    separate_rows(protein_id, sep = ",")

mRNA <- read.table(mRNA_file, sep = ",", col.names = c("scaffold", "start", "end", "comment")) %>%
    extract(comment, into = "product", regex = "product=([^;]+)", remove = F) %>%
    extract(comment, into = "dbxref", regex = "Dbxref=([^;]+)", remove = F) %>%
    extract(comment, into = "note", regex = "Note=([^;]+)", remove = F) %>%
    extract(comment, into = "ID", regex = "ID=([^;]+)", remove = F) %>%
    extract(comment, into = "protein_id", regex = "protein_id=ncbi:([^;]+)") %>%
    left_join(proteinortho, by = "protein_id") %>%
    filter(num_species > min_species) %>%
    group_by(scaffold) %>%
    summarize(n_homologs = n())

complement_gc <- read.table(complement_gc_file, col.names = c("scaffold", "compl_size", "compl_GC"))

clade_df <- lapply(clade_gc_files, read.table, col.names = c("scaffold", "size", "GC")) %>%
    setNames(clades) %>%
    bind_rows(.id = "clade") %>%
    left_join(MCP_scaf, by = c("scaffold", "clade")) %>%
    group_by(scaffold) %>%
    mutate(longest = clade[which.max(size)]) %>%
    pivot_wider(id_cols = c(scaffold, longest), names_from = clade, values_from = c(size, GC, MCP))

data <- read.table(coverage_file, col.names = c("scaffold", "scaf_size", "scaf_cov", "r"), fill = T) %>%
    filter(!is.na(r)) %>%
    select(-r) %>%
    left_join(complement_gc, by = "scaffold") %>%
    left_join(clade_df, by = "scaffold") %>%
    left_join(mRNA, by = "scaffold") %>%
    replace_na(list(compl_size = 0)) %>%
    mutate(scaf_cov_norm = scaf_cov / median(scaf_cov)) %>%
    mutate(viral_size = scaf_size - compl_size, viral_fract = viral_size / scaf_size) %>%
    filter(viral_size > 0)
write.csv(data, file = output_file, row.names = F, na = "")
