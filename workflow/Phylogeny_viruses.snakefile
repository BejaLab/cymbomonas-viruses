
small = [ 'Gezel', 'TVS', 'X' ]
big = [ 'Mesomimi', 'Dwarf' ]

small_genes = [ 'MCP', 'mCP_penton', 'Yrec', 'A32_ATPase', 'Tlr6F', 'PGVV05', 'PGVV09', 'YSL1_23', 'pPolB' ]

rule phylogeny:
    input:
        expand("output/phylogeny_viruses/{gene}-small.svg", gene = small_genes)

rule phylo_select:
    input:
        "output/mcl_genes_60.tsv"
    output:
        "analysis/phylogeny_viruses/{gene}-{group}.tsv"
    params:
        groups = lambda w: '|'.join(small if w.group == 'small' else big),
        regex = lambda w: '%s|%s' % (w.gene, w.gene.replace('_', ' '))
    conda:
        "envs/tools.yaml"
    shell:
        "csvgrep -t -c Gene -r '^({params.regex})$' {input} | csvgrep -c clade -r '^({params.groups})$' | csvcut -c ID,Genome | csvformat -T -K1 > {output}"

rule phylo_extract:
    input:
        tsv = "analysis/phylogeny_viruses/{gene}-{group}.tsv",
        fas = "analysis/hhblits_db/All_proteins.faa",
        fai = "analysis/hhblits_db/All_proteins.faa.fai"
    output:
        "analysis/phylogeny_viruses/{gene}-{group}.faa"
    conda:
        "envs/tools.yaml"
    shell:
        "cut -f1 {input.tsv} | xargs seqkit faidx {input.fas} | seqkit seq -o {output}"

rule phylo_mafft:
    input:
        "analysis/phylogeny_viruses/{gene}-{group}.faa"
    output:
        "analysis/phylogeny_viruses/{gene}-{group}.mafft"
    threads:
        20
    conda:
        "envs/tools.yaml"
    shell:
        "mafft --localpair --maxiterate 1000 --thread {threads} {input} > {output}"

rule phylo_trimal:
    input:
        "analysis/phylogeny_viruses/{gene}-{group}.mafft"
    output:
        "analysis/phylogeny_viruses/{gene}-{group}.trimal"
    params:
        gt = 0.9
    conda:
        "envs/tools.yaml"
    shell:
        "trimal -in {input} -out {output} -gt {params.gt}"

rule phylo_iqtree:
    input:
        "analysis/phylogeny_viruses/{gene}-{group}.trimal"
    output:
        "analysis/phylogeny_viruses/{gene}-{group}.treefile"
    params:
        seed = 123,
        B = 1000,
        prefix = "analysis/phylogeny_viruses/{gene}-{group}"
    threads:
        4
    conda:
        "envs/tools.yaml"
    shell:
        "iqtree2 -s {input} --prefix {params.prefix} -redo --alrt 1000 -B 1000 --seed {params.seed} -T {threads}"

rule phylo_ggtree:
    input:
        tree = "analysis/phylogeny_viruses/{gene}-{group}.treefile",
        tsv  = "analysis/phylogeny_viruses/{gene}-{group}.tsv",
        virus_metadata = "metadata/viruses.tsv",
        segment_metadata = "metadata/viral_segments.tsv",
        colors = "metadata/subclade_colors.txt"
    output:
        image = "output/phylogeny_viruses/{gene}-{group}.svg",
        jtree = "output/phylogeny_viruses/{gene}-{group}.jtree"
    params:
        outgroup_rooting = False
    conda:
        "envs/r.yaml"
    script:
        "scripts/ggtree-viruses.R"
