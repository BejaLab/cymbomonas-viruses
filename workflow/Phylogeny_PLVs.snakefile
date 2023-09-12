
rule phylogeny_PLVs:
    input:
        "output/phylogeny/PLVs-MCP.svg"

rule hmmsearch_extract_PLV_MCP:
    input:
        fasta = "analysis/PLVs/{genome}.faa",
        fai = "analysis/PLVs/{genome}.faa.seqkit.fai",
        tblout = "analysis/PLVs/{genome}.custom.tblout"
    output:
        "analysis/phylogeny/PLVs/viruses/{genome}-MCP.faa"
    params:
        evalue = 1e-8,
        regex = "major capsid protein",
        minlen = 400
    conda:
        "envs/tools.yaml"
    shell:
        "grep -i {params.regex:q} {input.tblout} | awk -ve={params.evalue} '$5<e{{print$3}}' | sort -u | xargs -r seqkit faidx -f {input.fasta} | seqkit seq -gm {params.minlen} -o {output}"

rule cat_cluster_MCPs:
    input:
        expand("analysis/locate_viruses/long_MCPs/{cluster}-{{clade}}.faa", cluster = clusters)
    output:
        "analysis/phylogeny/PLVs/clades/{clade}-MCP.faa"
    shell:
        "cat {input} > {output}"

rule cluster_MCPs_first:
    input:
        "metadata/curated/MCPs.faa"
    output:
        "analysis/phylogeny/PLVs/clades/{clade}-MCP_first.faa"
    conda:
        "envs/tools.yaml"
    shell:
        "(seqkit grep -rp @{wildcards.clade} {input} || true) | seqkit head -n1 -o {output}"

rule search_IMGVR:
    input:
        query = "analysis/phylogeny/PLVs/clades/{clade}-MCP_first.faa",
        db = "databases/IMGVR/IMGVR_all_nucleotides_renamed.fna"
    output:
        "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.tblastn"
    params:
        evalue = 1e-10
    conda:
        "envs/tools.yaml"
    threads:
        20
    shell:
        "tblastn -num_threads {threads} -query {input.query} -db {input.db} -outfmt 6 -out {output} -evalue {params.evalue}"

rule extract_IMGVR:
    input:
        tblastn = "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.tblastn",
        db = "databases/IMGVR/IMGVR_all_nucleotides_renamed.fna",
        tsv = "metadata/IMG.tsv"
    output:
        "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.fna"
    params:
        minlen = 15000
    conda:
        "envs/tools.yaml"
    shell:
        "cut -f2 {input.tblastn} | paste -sd, | xargs blastdbcmd -db {input.db} -outfmt %t$'\\t'%s -entry | sort -u | seqkit tab2fx | seqkit seq -gm {params.minlen} | seqkit grep -vrf <(grep Restricted {input.tsv} | cut -f1) -o {output}"

rule blast_IMGVR:
    input:
        db = "analysis/phylogeny/PLVs/clades/{clade}-MCP.cdhit",
        pdb = "analysis/phylogeny/PLVs/clades/{clade}-MCP.cdhit.pdb",
        query = "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.combined.faa"
    output:
        "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.blast"
    params:
        evalue = 1e-10
    conda:
        "envs/vcontact2.yaml"
    shell:
        "diamond blastp --query {input.query} --db {input.db} --evalue {params.evalue} --outfmt 6 --out {output} --more-sensitive"

rule get_IMGVR_MCPs:
    input:
        blast = "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.blast",
        faa = "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.combined.faa",
        fai = "analysis/phylogeny/PLVs/IMGVR_blast/{clade}.combined.faa.fai"
    output:
        "analysis/phylogeny/PLVs/IMGVR/{clade}.faa"
    params:
        bitscore = lambda w: PLV_scores[w.clade],
        minlen = 400
    conda:
        "envs/tools.yaml"
    shell:
        "awk '$12>{params.bitscore}' {input.blast} | cut -f1 | sort -u | xargs -r seqkit faidx {input.faa} | seqkit seq -m {params.minlen} -o {output}"

rule cdhit_MCPs:
    input:
        "analysis/phylogeny/PLVs/{source}/{basename}.faa"
    output:
        cdhit = "analysis/phylogeny/PLVs/{source}/{basename}.cdhit",
        clstr = "analysis/phylogeny/PLVs/{source}/{basename}.cdhit.clstr"
    wildcard_constraints:
        source = 'IMGVR|algae|clades|outgroups'
    params:
        c = 0.9
    conda:
        "envs/tools.yaml"
    shell:
        "cdhit -i {input} -o {output.cdhit} -c {params.c} -d 0"

rule blast_MCPs:
    input:
        db = "analysis/phylogeny/PLVs/clades/{clade}-MCP.cdhit",
        query = "analysis/algae/{genome}.faa"
    output:
        "analysis/phylogeny/PLVs/algae/{clade}-{genome}.blast"
    params:
        evalue = 1e-10
    conda:
        "envs/vcontact2.yaml"
    shell:
        "diamond blastp --query {input.query} --db {input.db} --evalue {params.evalue} --outfmt 6 --out {output} --more-sensitive"

rule get_algal_MCPs:
    input:
        blast = "analysis/phylogeny/PLVs/algae/{clade}-{genome}.blast",
        faa = "analysis/algae/{genome}.faa",
        fai = "analysis/algae/{genome}.faa.fai"
    output:
        "analysis/phylogeny/PLVs/algae/{clade}-{genome}.faa"
    params:
        bitscore = lambda w: PLV_scores[w.clade],
        minlen = 300
    conda:
        "envs/tools.yaml"
    shell:
        "awk '$12>{params.bitscore}' {input.blast} | cut -f1 | sort -u | xargs -r seqkit faidx {input.faa} | seqkit seq -m {params.minlen} -o {output}"

rule cat_PLV_MCPs_viruses:
    input:
        lambda w: expand("analysis/phylogeny/PLVs/viruses/{genome}-MCP.faa", genome = [ virus for x in PLV_outgroups for virus in virus_clades[x] if not virus.startswith('jcf') ])
    output:
        "analysis/phylogeny/PLVs/outgroups/MCP.faa"
    shell:
        "cat {input} > {output}"

rule mafft_MCPs:
    input:
        ingroup = lambda w: expand("analysis/phylogeny/PLVs/clades/{clade}-MCP.cdhit", clade = PLVs),
        outgroup = "analysis/phylogeny/PLVs/outgroups/MCP.cdhit",
        algae = lambda w: expand("analysis/phylogeny/PLVs/algae/{clade}-{genome}.cdhit", clade = PLVs, genome = algae),
        IMGVR = lambda w: expand("analysis/phylogeny/PLVs/IMGVR/{clade}.cdhit", clade = PLVs)
    output:
        "analysis/phylogeny/PLVs/PLVs.mafft"
    params:
        mafft = "--auto"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit rmdup -s {input} | mafft {params.mafft} --reorder - > {output}"

rule trimal_MCPs:
    input:
        "analysis/phylogeny/PLVs/PLVs.mafft"
    output:
        "analysis/phylogeny/PLVs/PLVs.trimal"
    params:
        gt = 0.9
    conda:
        "envs/tools.yaml"
    shell:
        "trimal -in {input} -out {output} -gt {params.gt}"

rule MCP_iqtree:
    input:
        "analysis/phylogeny/PLVs/PLVs.trimal"
    output:
        "analysis/phylogeny/PLVs/PLVs.treefile"
    params:
        seed = 123,
        B = 1000,
        prefix = "analysis/phylogeny/PLVs/PLVs"
    threads:
        4
    conda:
        "envs/tools.yaml"
    shell:
        "iqtree2 -s {input} --prefix {params.prefix} -redo --alrt 1000 -B {params.B} --seed {params.seed} -T {threads}"

rule faa_txt_viruses:
    input:
        "analysis/phylogeny/PLVs/viruses/{genome}-MCP.faa"
    output:
        "analysis/phylogeny/PLVs/viruses/{genome}-MCP.faa.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -ni {input} | awk -vg={wildcards.genome:q} '{{print g,$1}}' OFS='\\t' > {output}"

rule faa_txt_algae:
    input:
        "analysis/phylogeny/PLVs/algae/{clade}-{genome}.faa"
    output:
        "analysis/phylogeny/PLVs/algae/{clade}-{genome}.faa.txt"
    params:
        organism = lambda w: algae[w.genome]
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -ni {input} | awk -vg={params.organism:q} '{{print g,$1}}' OFS='\\t' > {output}"

rule PLV_phylo_ggtree:
    input:
        tree = "analysis/phylogeny/PLVs/PLVs.treefile",
        proteins = lambda w:
            expand("analysis/phylogeny/PLVs/viruses/{genome}-MCP.faa.txt", genome = [ virus for x in PLV_outgroups for virus in virus_clades[x]]) +
            expand("analysis/phylogeny/PLVs/algae/{clade}-{genome}.faa.txt", genome = algae, clade = PLVs),
        viruses  = "metadata/viruses.tsv",
        img = "metadata/IMG.tsv",
        clstr = lambda w: expand("analysis/phylogeny/PLVs/clades/{clade}-MCP.cdhit.clstr", clade = PLVs)
    output:
        image = "output/phylogeny/PLVs-MCP.svg",
        jtree = "output/phylogeny/PLVs-MCP.jtree"
    params:
        width = 4
    conda:
        "envs/r.yaml"
    script:
        "scripts/ggtree-PLVs.R"
