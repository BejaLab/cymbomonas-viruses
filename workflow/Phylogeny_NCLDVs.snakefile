
def NCLDV_segment_markers(wildcards):
    fasta_file = checkpoints.resolve_bed.get(cluster = clusters[0], clade = "NCLDV").output['fna']
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
        "grep -v '^#' {input.txt} | awk '$8>{params.score}' | cut -f1 -d' ' | head -n1 | xargs -r seqkit faidx {input.faa} | seqkit replace -p $ -r ' {wildcards.marker} {wildcards.genome}' -o {output} || true"

rule NCLDV_fetch:
    input:
        fna = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-NCLDV.fna",
        gff = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-NCLDV.combined.gff"
    output:
        "analysis/phylogeny/NCLDVs/faa/segments/{segment}.faa"
    conda:
        "envs/tools.yaml"
    shell:
        "grep -w {wildcards.segment:q} {input.gff} | gffread -g {input.fna} -w - | seqkit translate --trim | seqkit replace -p TCONS -r {wildcards.segment:q} -o {output}"

rule NCLDV_markers_mafft:
    input:
        expand("analysis/phylogeny/NCLDVs/hmmsearch/viruses/{genome}-{{marker}}.faa", genome = [ virus for x in NCLDVs for virus in virus_clades[x]]),
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
        nex = "analysis/phylogeny/NCLDVs/markers/NCLDV.nex",
        fasta = expand("analysis/phylogeny/NCLDVs/markers/{marker}.trimal", marker = NCLDV_markers)
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
        "iqtree2 -spp {input.nex} --prefix {params.prefix} -redo --alrt 1000 -B {params.B} --seed {params.seed} -T {threads}"

rule NCLDV_txt_viruses:
    input:
        "analysis/phylogeny/NCLDVs/faa/{source}/{genome}.faa"
    output:
        "analysis/phylogeny/NCLDVs/faa/{source}/{genome}.faa.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -ni {input} | awk -vg={wildcards.genome:q} '{{print g,$1}}' OFS='\\t' > {output}"

rule NCLDV_phylo_ggtree:
    input:
        tree = "analysis/phylogeny/NCLDVs/markers/NCLDV.treefile",
        proteins = lambda w: expand("analysis/phylogeny/NCLDVs/faa/{source}/{genome}.faa.txt", genome = [ virus for x in NCLDVs for virus in virus_clades[x]]),
        viruses  = "metadata/viruses.tsv",
        img = "metadata/IMG.tsv",
        clstr = "analysis/blank.txt"
    output:
        image = "output/phylogeny/NCLDVs.svg",
        jtree = "output/phylogeny/NCLDVs.jtree"
    conda:
        "envs/r.yaml"
    script:
        "scripts/ggtree-viruses.R"
