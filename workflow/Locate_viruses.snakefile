
from os.path import splitext
from Bio import SeqIO

rule locate_viruses:
    input:
        expand("analysis/locate_viruses/all_MCPs/{cluster}.bed", cluster = clusters)

rule MCPs_link:
    input:
        "metadata/curated/MCPs.faa"
    output:
        "analysis/locate_viruses/blastx/MCPs.faa"
    shell:
        "ln -rs {input} {output}"

rule MCPs_makeblastdb:
    input:
        "analysis/locate_viruses/blastx/MCPs.faa"
    output:
        "analysis/locate_viruses/blastx/MCPs.faa.pdb"
    conda:
        "envs/tools.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

rule blast_complete:
    input:
        db = "analysis/locate_viruses/blastx/MCPs.faa",
        pdb = "analysis/locate_viruses/blastx/MCPs.faa.pdb",
        fasta = "genomes/viruses.fna"
    output:
        "analysis/locate_viruses/blastx/viruses.blastx"
    params:
        evalue = 1e-100
    threads:
        10
    conda:
        "envs/tools.yaml"
    shell:
        "blastx -num_threads {threads} -db {input.db} -query {input.fasta} -evalue {params.evalue} -max_target_seqs 1 -outfmt 6 -out {output}"

rule get_clade:
    input:
        fasta = "genomes/viruses.fna",
        blastx = "analysis/locate_viruses/blastx/viruses.blastx"
    output:
        "analysis/locate_viruses/clades/{clade}.fna"
    shell:
        "grep @{wildcards.clade} {input.blastx} | cut -f1 | seqkit grep -f- -o {output} {input.fasta}"

rule clusters_link:
    input:
        "clusters/{cluster}.fasta"
    output:
        "analysis/locate_viruses/clusters/{cluster}.fna"
    shell:
        "ln -rs {input} {output}"

rule blastn:
    input:
        query = "analysis/locate_viruses/clades/{clade}.fna",
        db  = "analysis/locate_viruses/clusters/{cluster}.fna",
        ndb = "analysis/locate_viruses/clusters/{cluster}.fna.ndb"
    output:
        "analysis/locate_viruses/blastn/{cluster}-{clade}.xml"
    params:
        evalue = 1e-10
    conda:
        "envs/tools.yaml"
    shell:
        "blastn -db {input.db} -query {input.query} -evalue {params.evalue} -outfmt 5 -out {output}"

rule parse_blastn:
    input:
        xml = "analysis/locate_viruses/blastn/{cluster}-{clade}.xml",
        fasta = "analysis/locate_viruses/clusters/{cluster}.fna",
        query = "analysis/locate_viruses/clades/{clade}.fna"
    output:
        "analysis/locate_viruses/blastn/{cluster}-{clade}.gff"
    params:
        query_gap = 1000,
        sbjct_gap_1 = 10000,
        sbjct_gap_2 = 1000,
        chain_len = 1000,
        min_neighbor_len = 500
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/blastn.py"

rule simplify_gff:
    input:
        "analysis/locate_viruses/blastn/{cluster}-{clade}.gff"
    output:
        "analysis/locate_viruses/blastn/{cluster}-{clade}.gff3"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/simplify_gff.py"

rule flatten_gff:
    input:
        "analysis/locate_viruses/blastn/{cluster}-{clade}.gff3"
    output:
        "analysis/locate_viruses/segments/{cluster}-{clade}.gff3"
    params:
        overlap = 10
    conda:
        "envs/tools.yaml"
    shell:
        """
        sort -k1,1 -k4,4n {input} | bedtools merge -s -d -{params.overlap} -c 7 -o distinct | awk '{{print$1,"blastn","fragment",$2+1,$3,".",$4,".",sprintf("ID=frag-%s;", NR)}}' OFS=\\\\t > {output}
        """

rule gff3_to_fna:
    input:
        gff = "analysis/locate_viruses/segments/{cluster}-{clade}.gff3",
        fna = "analysis/locate_viruses/clusters/{cluster}.fna"
    output:
        "analysis/locate_viruses/segments/{cluster}-{clade}.fna"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/gff_fna.py"

rule segment_hmmsearch:
    input:
        hmm = "analysis/hmm/custom.MCPs.hmm",
        faa = "analysis/locate_viruses/segments/{cluster}-{clade}.combined.faa"
    output:
        "analysis/locate_viruses/MCP/{cluster}-{clade}.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmsearch --tblout {output} -o /dev/null {input.hmm} {input.faa}"

rule segment_hmmsearch_fasta:
    input:
        faa = "analysis/locate_viruses/segments/{cluster}-{clade}.combined.faa",
        fai = "analysis/locate_viruses/segments/{cluster}-{clade}.combined.faa.fai",
        txt = "analysis/locate_viruses/segments/{cluster}-{clade}.combined.custom.tblout"
    output:
        "analysis/locate_viruses/MCP/{cluster}-{clade}.faa"
    params:
        evalue = 1e-5,
        regex = "major capsid protein"
    conda:
        "envs/tools.yaml"
    shell:
        "grep -i {params.regex:q} {input.txt} | awk '$5<{params.evalue}{{print$3}}' | xargs seqkit faidx {input.faa} | seqkit replace -p TCONS -r {wildcards.cluster}_{wildcards.clade} > {output}"

rule segment_hmmsearch_fasta_cat:
    input:
        "analysis/locate_viruses/MCP/{cluster}-{clade}.faa"
    output:
        "analysis/locate_viruses/long_MCPs/{cluster}-{clade}.faa"
    params:
        m = 400
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -gm{params.m} {input} | seqkit rmdup -so {output}"

rule segment_hmmsearch_fasta_cdhit:
    input:
        "analysis/locate_viruses/long_MCPs/{cluster}-{clade}.faa"
    output:
        "analysis/locate_viruses/long_MCPs/{cluster}-{clade}.cdhit"
    params:
        c = 0.95
    conda:
        "envs/tools.yaml"
    shell:
        "cdhit -i {input} -o {output} -d 0 -c {params.c}"

rule MCP_blast:
    input:
        query = "analysis/locate_viruses/long_MCPs/{cluster}-{clade}.cdhit",
        fna = "analysis/locate_viruses/clusters/{cluster}.fna",
        ndb = "analysis/locate_viruses/clusters/{cluster}.fna.ndb"
    output:
        "analysis/locate_viruses/long_MCPs/{cluster}-{clade}.tblastn.outfmt6"
    params:
        evalue = 1e-10
    conda:
        "envs/tools.yaml"
    threads:
        20
    shell:
        "tblastn -num_threads {threads} -query {input.query} -db {input.fna} -evalue {params.evalue} -max_target_seqs 1000000 -outfmt 6 -out {output}"

rule outfmt6_to_bed:
    input:
        "analysis/locate_viruses/long_MCPs/{cluster}-{clade}.tblastn.outfmt6"
    output:
        "analysis/locate_viruses/long_MCPs/{cluster}-{clade}.tblastn.outfmt6.bed"
    conda:
        "envs/tools.yaml"
    shell:
        """
        awk -v c={wildcards.clade} '{{n=$2;b=$9;e=$10;s="+"}}b>e{{s="-";b=$10;e=$9}}{{print n,b-1,e,$3,c,s}}' OFS=\\\\t {input} | sort -k1,1 -k2,2n | bedtools merge -s -c 4,5,6 -o max,distinct,distinct > {output}
        """

rule merge_outfmt6_to_bed:
    input:
        expand("analysis/locate_viruses/long_MCPs/{{cluster}}-{clade}.tblastn.outfmt6.bed", clade = clades)
    output:
        "analysis/locate_viruses/all_MCPs/{cluster}.bed"
    conda:
        "envs/tools.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} | bedtools merge -s -c 4,5,6 -o collapse,collapse,distinct | awk '{{split($4,s,",");split($5,c,",");S=0;for(i in s) if(s[i]>S){{S=s[i];C=c[i]}};print$1,$2,$3,S,C,$6}}' OFS=\\\\t > {output}
        """
