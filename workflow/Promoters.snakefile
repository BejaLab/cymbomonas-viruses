
promoter_re = "[ATW][ATW][ATW][ATW][ATW]TG[ATW]"

rule promoters:
    input:
        "output/fimo/Gezel-14T.gff",
        expand("output/fimo/{segment}.gff", segment = Gezel_segments),
        "output/promoters.svg"

rule flanks:
    input:
        "analysis/PLVs/{genome}.gbk"
    output:
        "analysis/promoters/{genome}.fna"
    params:
        min_len = 100,
        max_len = 150
    conda:
        "envs/bioperl.yaml"
    shell:
        "perl workflow/helpers/biotags.pl -i {input} -p CDS -t seq-{params.max_len} | awk '$1' | nl -w1 | seqkit tab2fx | seqkit seq -gm{params.min_len} -o {output}"

rule meme:
    input:
        "analysis/promoters/{genome}.fna"
    output:
        "analysis/promoters/{genome}/meme/meme.xml",
        "analysis/promoters/{genome}/meme/meme.txt"
    params:
        outdir = "analysis/promoters/{genome}/meme/",
        minw = 6,
        maxw = 16,
        nmotifs = 10,
        minsites = 10
    conda:
        "envs/meme.yaml"
    shell:
        "meme -minw {params.minw} -maxw {params.maxw} -nmotifs {params.nmotifs} -minsites {params.minsites} -maxsites Inf -dna -oc {params.outdir} {input}"

rule meme_extract:
    input:
        "analysis/promoters/{genome}/meme/meme.txt"
    output:
        "analysis/promoters/{genome}/meme/meme.extract.txt"
    params:
        motif_re = promoter_re
    shell:
        """
        awk -vm='{params.motif_re}' '$1=="MOTIF"{{s=$2!~m}}!s' {input} > {output}
        """

rule gbk2fna:
    input:
        "metadata/Gezel-14T.gbk"
    output:
        "analysis/Gezel-14T/Gezel-14T.fna"
    conda:
        "envs/emboss.yaml"
    shell:
        "seqret -supper1 -filter -outseq {output} {input}"

rule fimo:
    input:
        meme = "analysis/promoters/PgV-16T/meme/meme.extract.txt",
        fasta = "analysis/Gezel-14T/Gezel-14T.fna"
    output:
        "analysis/fimo/Gezel-14T/fimo.tsv"
    params:
        out_dir= "analysis/fimo/Gezel-14T"
    conda:
        "envs/meme.yaml"
    shell:
        "fimo --oc {params.out_dir} {input.meme} {input.fasta}"

rule fimo_segments:
    input:
        meme = "analysis/promoters/PgV-16T/meme/meme.extract.txt",
        fasta = "analysis/PLV_segments/{segment}.fna"
    output:
        "analysis/fimo/{segment}/fimo.tsv"
    params:
        out_dir= "analysis/fimo/{segment}"
    conda:
        "envs/meme.yaml"
    shell:
        "fimo --oc {params.out_dir} {input.meme} {input.fasta}"

rule fimo_filter:
    input:
        "analysis/fimo/{item}/fimo.tsv"
    output:
        "output/fimo/{item}.gff"
    params:
        q_thresh = 0.05,
        motif_re = promoter_re
    conda:
        "envs/r.yaml"
    script:
        "scripts/fimo.R"

rule plot_meme:
    input:
        "analysis/promoters/{genome}/meme/meme.xml"
    output:
        "analysis/promoters/{genome}/meme/meme.svg"
    params:
        motif_res = [ promoter_re, "TCCGGA" ],
        meme_name = "{genome}",
        max_len   = 150
    conda:
        "envs/r.yaml"
    script:
        "scripts/meme.R"

rule merge_meme:
    input:
        expand("analysis/promoters/{virus}/meme/meme.svg", virus = Mesomimi_viruses)
    output:
        "output/promoters.svg"
    shell:
        "python workflow/scripts/svg_stack.py --direction=v {input} > {output}"

