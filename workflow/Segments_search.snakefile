
rule search_segments:
    input:
        expand("analysis/segments/{genome}/segments.gbk", genome = genomes)

rule getorf:
    input:
        "genomes/{genome}.fna"
    output:
        "analysis/getorf/{genome}.faa"
    conda:
        "envs/emboss.yaml"
    shell:
        "getorf -minsize 100 -maxsize 100000 -filter {input} > {output}"

rule hmmsearch:
    input:
        fasta = "analysis/getorf/{genome}.faa", hmm = "hmm/{hmm}.hmm"
    output:
        "analysis/hmmsearch/{genome}-{hmm}.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmsearch -o /dev/null --tblout {output} {input.hmm} {input.fasta}"

rule extract_hmmsearch:
    input:
        fasta = "analysis/getorf/{genome}.faa",
        fai = "analysis/getorf/{genome}.faa.seqkit.fai",
        txt = "analysis/hmmsearch/{genome}-{hmm}.txt"
    output:
        "analysis/hmmsearch/{genome}-{hmm}.faa"
    params:
        evalue = 0.001
    conda:
        "envs/tools.yaml"
    shell:
        """
        awk '!/^#/&&$5<{params.evalue}{{print $1,$3}}' OFS=\\\\t {input.txt} | \
            parallel --colsep \\\\t seqkit faidx -f {input.fasta} {{1}} \| seqkit replace -p "'$'" -r '" {{2}}"' | seqkit  seqkit seq -gm 100
        """

rule parse_hmmsearch:
    input:
        expand("analysis/hmmsearch/{{genome}}-{hmm}.txt", hmm = hmms)
    output:
        # figure = "images/segments_{genome}.pdf",
        bed = "analysis/segments/{genome}/segments.bed"
    params:
        e_value_threshold = 0.001,
        distance_threshold = 100000,
        genes_threshold = 2
    conda:
        "envs/r.yaml"
    script:
        "scripts/parse_hmmsearch.R"

rule segments_extract:
    input:
        genome = "genomes/{genome}.fna", bed = "analysis/segments/{genome}/segments.bed"
    output:
        "analysis/segments/{genome}/segments.fna"
    params:
        flanks = 50000
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit subseq --up-stream {params.flanks} --down-stream {params.flanks} --bed {input.bed} {input.genome} | cut -f1 -d: > {output}"

checkpoint segments_split:
    input:
        "analysis/segments/{genome}/segments.fna"
    output:
        directory("analysis/segments/{genome}/segments_split")
    conda:
        "envs/tools.yaml"
    shell:
        """
        seqkit split -iO {output} {input}
        rename s/segments.id_// {output}/*.fna
        """

rule segments_genemarks:
    input:
        "analysis/segments/{genome}/segments_split/{segment}.fna"
    output:
        "analysis/segments/{genome}/segments_split/{segment}.genemarks.gff"
    shadow:
        "minimal"
    shell:
        "gmsn.pl --format GFF3 --virus --output {output} {input}"

rule segments_prodigal:
    input:
        "analysis/segments/{genome}/segments_split/{segment}.fna"
    output:
        "analysis/segments/{genome}/segments_split/{segment}.prodigal.gff"
    shadow:
        "minimal"
    conda:
        "envs/prokka.yaml"
    shell:
        "prodigal -i {input} -m -g 1 -p meta -f gff > {output}"

rule segments_genemarks_cat:
    input:
        "analysis/segments/{genome}/segments_split/{segment}.genemarks.gff",
        "analysis/segments/{genome}/segments_split/{segment}.prodigal.gff"
    output:
        "analysis/segments/{genome}/segments_split/{segment}.combined.gtf"
    params:
        outprefix = "analysis/segments/{genome}/segments_split/{segment}",
        cprefix = lambda wildcards:
            wildcards.genome.replace('-','') + '_' + sha256(wildcards.segment.encode('utf-8')).hexdigest()[0:6]
    conda:
        "envs/tools.yaml" # NB: conda gives 0.11.2
    shell:
        "gffcompare -p {params.cprefix} -CTo {params.outprefix} {input}"

rule segments_gffread:
    input:
        gtf = "analysis/segments/{genome}/segments_split/{segment}.combined.gtf",
        fna = "analysis/segments/{genome}/segments_split/{segment}.fna",
        fai = "analysis/segments/{genome}/segments_split/{segment}.fna.fai"
    output:
        gff = "analysis/segments/{genome}/segments_split/{segment}.gff",
        faa = "analysis/segments/{genome}/segments_split/{segment}.faa"
    conda:
        "envs/tools.yaml"
    shell:
        "gffread -g {input.fna} -w - -o {output.gff} {input.gtf} | seqkit translate --trim -o {output.faa}"

def aggregate_segments(wildcards, suffix):
    split_dir = checkpoints.segments_split.get(**wildcards).output[0]
    fna = path.join(split_dir, "{segment}.fna")
    segments ,= glob_wildcards(fna)
    files = path.join(split_dir, suffix)
    return expand(files, segment = segments)

rule segments_gff_cat:
    input:
        lambda wildcards: aggregate_segments(wildcards, "{segment}.gff")
    output:
        "analysis/segments/{genome}/segments.gff"
    shell:
        "cat {input} > {output}"

rule segments_faa_cat:
    input:
        lambda wildcards: aggregate_segments(wildcards, "{segment}.faa")
    output:
        "analysis/segments/{genome}/segments.faa"
    shell:
        "cat {input} > {output}"

rule segments_gvog:
    input:
        hmm = "databases/GVDB/output/GVOG.hmm",
        fasta = "analysis/segments/{genome}/segments.faa"
    output:
        tblout = "analysis/segments/{genome}/segments.gvog.tblout",
        domtblout = "analysis/segments/{genome}/segments.gvog.domtblout"
    threads:
        4
    conda:
        "envs/tools.yaml"
    shell:
        "hmmscan --cpu {threads} -o /dev/null --tblout {output.tblout} --domtblout {output.domtblout} {input.hmm} {input.fasta}"

rule segments_bellas:
    input:
        hmm = "databases/Bellas_Sommaruga/output/PCs.hmmdb",
        fasta = "analysis/segments/{genome}/segments.faa"
    output:
        tblout = "analysis/segments/{genome}/segments.bellas.tblout",
        domtblout = "analysis/segments/{genome}/segments.bellas.domtblout"
    threads:
        4
    conda:
        "envs/tools.yaml"
    shell:
        "hmmscan --cpu {threads} -o /dev/null --tblout {output.tblout} --domtblout {output.domtblout} {input.hmm} {input.fasta}"

rule segments_vogdb:
    input:
        hmm = "databases/vogdb/output/vog.hmmdb",
        fasta = "analysis/segments/{genome}/segments.faa"
    output:
        tblout = "analysis/segments/{genome}/segments.vogdb.tblout",
        domtblout = "analysis/segments/{genome}/segments.vogdb.domtblout"
    threads:
        4
    conda:
        "envs/tools.yaml"
    shell:
        "hmmscan --cpu {threads} -o /dev/null --tblout {output.tblout} --domtblout {output.domtblout} {input.hmm} {input.fasta}"

rule segments_hmmsearch:
    input:
        fasta = "analysis/segments/{genome}/segments.faa", hmm = "hmm/{hmm}.hmm"
    output:
        "analysis/segments/{genome}/segments-hmm-{hmm}.tblout"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmsearch -o /dev/null --tblout {output} {input.hmm} {input.fasta}"

rule segments_pfam:
    input:
        "analysis/segments/{genome}/segments.faa"
    output:
        "analysis/segments/{genome}/segments.pfam.txt"
    params:
        pfam_dir = "databases/Pfam"
    threads:
        4
    conda:
        "envs/bioperl.yaml"
    shell:
        "pfam_scan.pl -fasta {input} -dir {params.pfam_dir} -outfile {output} -cpu {threads}"

rule segments_blast:
    input:
        fna = "analysis/segments/{genome}/segments.fna",
        ndb = "analysis/segments/{genome}/segments.fna.ndb"
    output:
        "analysis/segments/{genome}/segments.fna.blastn"
    conda:
        "envs/tools.yaml"
    shell:
        "blastn -query {input.fna} -db {input.fna} -outfmt 6 -out {output} -evalue 1e-20"

rule segments_genbank:
    input:
        fna = "analysis/segments/{genome}/segments.fna",
        gff = "analysis/segments/{genome}/segments.gff",
        hmm_dbs = [ "databases/Pfam/Pfam-A.hmm" ] + expand("hmm/{hmm}.hmm", hmm = hmms),
        blastn = "analysis/segments/{genome}/segments.fna.blastn",
        hmmscan = [
            "analysis/segments/{genome}/segments.gvog.tblout",
            "analysis/segments/{genome}/segments.vogdb.tblout",
            "analysis/segments/{genome}/segments.bellas.tblout"
        ],
        markers = expand("analysis/segments/{{genome}}/segments-hmm-{hmm}.tblout", hmm = hmms),
        pfam = "analysis/segments/{genome}/segments.pfam.txt"
    output:
        "analysis/segments/{genome}/segments.gbk"
    params:
        coding_feature = "exon",
        evalue = 1e-5,
        min_repeat_len = 70,
        min_repeat_gap = 1000,
        min_repeat_id  = 90
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/gff2gbk.py"
