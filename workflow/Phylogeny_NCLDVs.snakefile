
NCLDV_genomes = [ virus for x in NCLDVs for virus in virus_clades[x] if viruses[virus]['representative'] != 'y' ]

def NCLDV_segment_markers(wildcards):
    fasta_file = checkpoints.gff3_to_fna.get(cluster = clusters[0], clade = "NCLDV").output[0]
    recs = SeqIO.parse(fasta_file, "fasta")
    min_len = 80000
    segments = [ rec.id for rec in recs if len(rec.seq) > min_len ]
    return expand("analysis/phylogeny/NCLDVs/hmmsearch/segments/{segment}-{marker}.faa", segment = segments, marker = wildcards.marker)

rule phylogeny_NCLDVs:
    input:
        "output/phylogeny/NCLDVs.svg", "output/phylogeny/NCLDVs.jtree"

rule fetch_ncldv_hmm:
    input:
        "databases/GVDB/output/GVOG.hmm"
    output:
        "analysis/phylogeny/NCLDVs/hmm/{marker}.hmm"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmfetch {input} {wildcards.marker} > {output}"

rule ncldv_link_faa:
    input:
        "analysis/PLVs/{genome}.faa"
    output:
        "analysis/phylogeny/NCLDVs/faa/viruses/{genome}.faa"
    shell:
        "ln -sr {input} {output}"

rule ncldv_hmmsearch:
    input:
        hmm = "analysis/phylogeny/NCLDVs/hmm/{marker}.hmm",
        faa = "analysis/phylogeny/NCLDVs/faa/{source}/{genome}.faa"
    output:
        "analysis/phylogeny/NCLDVs/hmmsearch/{source}/{genome}-{marker}.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmsearch -o /dev/null --domtblout {output} {input.hmm} {input.faa}"

rule ncldv_hmmsearch_extract:
    input:
        faa = "analysis/phylogeny/NCLDVs/faa/{source}/{genome}.faa",
        fai = "analysis/phylogeny/NCLDVs/faa/{source}/{genome}.faa.fai",
        txt = "analysis/phylogeny/NCLDVs/hmmsearch/{source}/{genome}-{marker}.txt"
    output:
        "analysis/phylogeny/NCLDVs/hmmsearch/{source}/{genome}-{marker}.faa"
    params:
        score = lambda w: NCLDV_markers[w.marker]
    conda:
        "envs/tools.yaml"
    shell:
        "grep -v '^#' {input.txt} | awk '$8>{params.score}' | cut -f1 -d' ' | head -n1 | xargs -r seqkit faidx {input.faa} | seqkit replace -p ^ -r '{wildcards.genome} ' -o {output} || true"

rule NCLDV_fetch:
    input:
        fna = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-NCLDV.fna",
        gff = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-NCLDV.combined.gff"
    output:
        "analysis/phylogeny/NCLDVs/faa/segments/{segment}.faa"
    conda:
        "envs/tools.yaml"
    shell:
        "grep -w {wildcards.segment:q} {input.gff} | gffread -g {input.fna} -w - | seqkit translate --trim -o {output}"

rule NCLDV_markers_mafft:
    input:
        expand("analysis/phylogeny/NCLDVs/hmmsearch/viruses/{genome}-{{marker}}.faa", genome = NCLDV_genomes),
        NCLDV_segment_markers
    output:
        "analysis/phylogeny/NCLDVs/markers/{marker}.mafft"
    conda:
        "envs/tools.yaml"
    shell:
        "cat {input} | mafft --auto - > {output}"

rule NCLDV_markers_trimal:
    input:
        "analysis/phylogeny/NCLDVs/markers/{marker}.mafft"
    output:
        "analysis/phylogeny/NCLDVs/markers/{marker}.trimal"
    params:
        gt = 0.9
    conda:
        "envs/tools.yaml"
    shell:
        "trimal -gt {params.gt} -in {input} -out {output}"

rule NCLDV_markers_nexus:
    input:
        expand("analysis/phylogeny/NCLDVs/markers/{marker}.trimal", marker = NCLDV_markers)
    output:
        "analysis/phylogeny/NCLDVs/markers/NCLDV.nex"
    conda:
        "envs/tools.yaml"
    shell:
        """
        (
            echo "#nexus"
            echo "begin sets;"
            parallel echo '"charset marker{{#}} = {{}}:*;"' ::: {input:q}
            echo "end;"
        ) > {output}
        """
 
rule NCLDV_iqtree:
    input:
        "analysis/phylogeny/NCLDVs/markers/NCLDV.nex",
    output:
        "analysis/phylogeny/NCLDVs/markers/NCLDV.treefile"
    params:
        seed = 123,
        B = 1000,
        prefix = "analysis/phylogeny/NCLDVs/markers/NCLDV"
    threads:
        4
    conda:
        "envs/tools.yaml"
    shell:
        "iqtree2 -spp {input} --prefix {params.prefix} -redo --alrt 1000 -B {params.B} --seed {params.seed} -T {threads}"

rule NCLDV_txt_viruses:
    input:
        expand("analysis/phylogeny/NCLDVs/markers/{marker}.mafft", marker = NCLDV_markers)
    output:
        "analysis/phylogeny/NCLDVs/markers/names.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -ni {input} | sort -u | awk '{{print $1,$1}}' OFS='\\t' > {output}"

rule NCLDV_phylo_ggtree:
    input:
        tree = "analysis/phylogeny/NCLDVs/markers/NCLDV.treefile.madroot",
        proteins = "analysis/phylogeny/NCLDVs/markers/names.txt",
        viruses  = "metadata/viruses.tsv",
        markers  = "metadata/NCLDV_markers.txt",
        img = "metadata/IMG.tsv",
        clstr = "analysis/blank.txt",
        lens =
            expand("analysis/PLVs/{genome}.lens.txt", genome = NCLDV_genomes) +
            expand("analysis/locate_viruses/segments/{cluster}-NCLDV.lens.txt", cluster = clusters),
        alns = expand("analysis/phylogeny/NCLDVs/markers/{marker}.mafft", marker = NCLDV_markers)
    output:
        image = "output/phylogeny/NCLDVs.svg",
        jtree = "output/phylogeny/NCLDVs.jtree"
    params:
        width = 7
    conda:
        "envs/r.yaml"
    script:
        "scripts/ggtree-viruses.R"
