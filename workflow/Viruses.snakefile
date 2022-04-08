
from pandas import read_csv, isnull

segments_file = 'metadata/viral_elements.tsv'
viruses_file  = 'metadata/viruses.tsv'
bellas_file   = 'metadata/viruses_bellas.tsv'
yutin_file    = 'metadata/viruses_yutin.tsv'

# search = [ "autoblast", "blastp", "diamond", "usearch", "ublast", "lastp", "rapsearch", "topaz", "blatp", "mmseqsp" ]
search = [ "diamond" ]

segments = {}
data = read_csv(segments_file, sep = '\t')
for i, row in data.iterrows():
    if isnull(row['fragment']):
        short = row['short']
        segments[short] = row

viruses = {}
data = read_csv(viruses_file, sep = '\t')
for i, row in data.iterrows():
    if isnull(row['fragment']):
        short = row['short']
        viruses[short] = row

bellas = {}
data = read_csv(bellas_file, sep = '\t')
for i, row in data.iterrows():
    if isnull(row['fragment']):
        short = row['short']
        bellas[short] = row

yutin = {}
data = read_csv(yutin_file, sep = '\t')
for i, row in data.iterrows():
    if isnull(row['fragment']):
        short = row['short']
        yutin[short] = row

rule all:
    input:
        expand("analysis/PLVs/orthogroups/{search}", search = search)

rule ref_element:
    output:
        "analysis/PLVs/viruses/{segment}.gbk"
    params:
        id = lambda w: viruses[w.segment]['accession']
    resources:
        ncbi = 1
    shell:
        "efetch -db nuccore -id {params.id} -format gb > {output}"

rule extract_element:
    input:
        lambda w: "analysis/segments/%s/segments.gbk" % segments[w.segment]['genome']
    output:
        "analysis/PLVs/elements/{segment}.gbk"
    params:
        flank = 50000,
        segment = lambda w: segments[w.segment]
    script:
        "scripts/extract_element.py"

# Workaround for BQ2
rule extract_faa_BQ2:
    input:
        "analysis/PLVs/viruses/CpV-BQ2.gbk"
    output:
        "analysis/PLVs/viruses/CpV-BQ2.faa"
    shell:
        "sed s/misc_feature/CDS/ {input} | biotags.pl -i - -f genbank -p CDS -t note,seq_translate | awk -F'[ ]' '{{print$NF}}' | seqkit tab2fx -o {output}"

rule extract_faa_bellas:
    input:
        "databases/Bellas_Sommaruga/input/All_proteins.faa"
    output:
        "analysis/PLVs/bellas/{segment}.faa"
    params:
        search = lambda w: bellas[w.segment]['search']
    shell:
        "seqkit grep -rp {params.search} -o {output} {input}"

rule extract_gbk_yutin:
    input:
        "databases/Yutin2015.gbk"
    output:
        "analysis/PLVs/yutin/{segment}.gbk"
    shell:
        "awk -vl={wildcards.segment} '/^LOCUS/{p=$2==l}p' {input} | sed -E 's/\\x93|\\x94/\"/g'"

rule extract_faa_yutin:
    input:
        "analysis/PLVs/yutin/{segment}.gbk"
    output:
        "analysis/PLVs/yutin/{segment}.faa"
    shell:
        "biotags.pl -i {input} -p CDS -t protein_id,translation | awk '$1' | seqkit tab2fx -o {output}"

rule extract_faa_viruses:
    input:
        "analysis/PLVs/viruses/{segment}.gbk"
    output:
        "analysis/PLVs/viruses/{segment}.faa"
    shell:
        "biotags.pl -i {input} -p CDS -t protein_id,translation | awk '$1' | seqkit tab2fx -o {output}"

rule extract_faa_elements:
    input:
        "analysis/PLVs/elements/{segment}.gbk"
    output:
        "analysis/PLVs/elements/{segment}.faa"
    shell:
        "biotags.pl -i {input} -p CDS -t locus_tag,translation | awk '$1' | seqkit tab2fx -o {output}"

rule proteinortho:
    input:
        expand("analysis/PLVs/elements/{segment}.faa", segment = segments.keys()),
        expand("analysis/PLVs/viruses/{segment}.faa",  segment = viruses.keys()),
        expand("analysis/PLVs/bellas/{segment}.faa",  segment = bellas.keys())
    output:
        "analysis/PLVs/proteinortho/{search}/proteinortho.proteinortho.tsv"
    shadow:
        "shallow"
    params:
        evalue = 1e-5,
        conn   = 0.1,
        outdir = "analysis/PLVs/proteinortho/{search}/"
    threads:
        workflow.cores
    shell:
        """
        proteinortho -project=proteinortho -p={wildcards.search} -cpus={threads} -conn={params.conn} -e={params.evalue} {input}
        mv proteinortho.* {params.outdir}
        """

checkpoint proteinortho_collect:
    input:
        "analysis/PLVs/proteinortho/{search}/proteinortho.proteinortho.tsv",
        expand("analysis/PLVs/elements/{segment}.faa", segment = segments.keys()),
        expand("analysis/PLVs/viruses/{segment}.faa",  segment = viruses.keys()),
        expand("analysis/PLVs/bellas/{segment}.faa",  segment = bellas.keys())
    output:
        directory("analysis/PLVs/orthogroups/{search}")
    shadow:
        "shallow"
    shell:
        "perl workflow/scripts/proteinortho_grab_proteins.pl -t -s -output_dir={output} -output_file=sequences.faa -exact {input}"
