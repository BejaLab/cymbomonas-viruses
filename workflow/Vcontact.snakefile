
rule vcontact:
    input:
        "analysis/vcontact2/results",
        "output/vcontact2_clustering.svg"

rule vcontact2_cat:
    input:
        expand("analysis/PLVs/{genome}.faa", genome = Gezel_viruses)
    output:
        "analysis/vcontact2/proteins.faa"
    shell:
        "cat {input} > {output}"

rule vcontact2_map:
    input:
        expand("analysis/PLVs/{genome}.faa", genome = Gezel_viruses)
    output:
        "analysis/vcontact2/proteins.csv"
    conda:
        "envs/tools.yaml"
    shell:
        """
        parallel --tagstring {{/.}} seqkit seq -ni {{}} ::: {input} | awk -vOFS=, 'BEGIN{{print"protein_id","contig_id","keywords"}}{{print$2,$1,""}}' > {output}
        """

rule vcontact2_run:
    input:
        fasta = "analysis/vcontact2/proteins.faa",
        mapping = "analysis/vcontact2/proteins.csv"
    output:
        outdir = directory("analysis/vcontact2/results"),
        network  = "analysis/vcontact2/results/c1.ntw",
        genomes  = "analysis/vcontact2/results/genome_by_genome_overview.csv",
        profiles = "analysis/vcontact2/results/vConTACT_profiles.csv",
    conda:
        "envs/vcontact2.yaml"
    params:
        db = "None",
        rel_mode = "Diamond",
        pcs_mode = "MCL",
        vcs_mode = "ClusterONE"
    threads:
        workflow.cores
    shell:
        "vcontact2 --threads {threads} --db {params.db} --raw-proteins {input.fasta} --rel-mode {params.rel_mode} --proteins-fp {input.mapping} --pcs-mode {params.pcs_mode} --vcs-mode {params.vcs_mode} --c1-bin $CONDA_PREFIX/lib/cluster_one-v1.0.jar --output-dir {output.outdir}"

rule vcontact2_plot:
    input:
        network  = "analysis/vcontact2/results/c1.ntw",
        genomes  = "analysis/vcontact2/results/genome_by_genome_overview.csv",
        profiles = "analysis/vcontact2/results/vConTACT_profiles.csv",
        segments = "metadata/viral_segments.tsv",
        viruses  = "metadata/viruses.tsv",
        colors   = "metadata/subclade_colors.txt"
    output:
        clustering = "output/vcontact2_clustering.svg"
    conda:
        "envs/r.yaml"
    script:
        "scripts/vcontact2.R"
