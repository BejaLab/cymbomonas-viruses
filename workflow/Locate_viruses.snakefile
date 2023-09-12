
from os.path import splitext
from Bio import SeqIO

rule locate_viruses:
    input:
        #expand("output/locate_viruses_{cluster}.svg", cluster = clusters),
        #expand("output/scaffold_{scaffold}.svg", scaffold = [ "jcf7180000139292", "jcf7180000174485" ]),
        expand("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.stringtie.faa", clade = clades),
        "output/total_coverage.pdf"

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
        query_gap   = lambda w: 1000  if w.clade != 'NCLDV' else 3000,
        sbjct_gap_1 = lambda w: 10000 if w.clade != 'NCLDV' else 30000,
        sbjct_gap_2 = lambda w: 1000  if w.clade != 'NCLDV' else 3000,
        chain_len   = lambda w: 1000  if w.clade != 'NCLDV' else 3000,
        min_neighbor_len = 500
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/blastn.py"

rule simplify_gff:
    input:
        "analysis/locate_viruses/blastn/{cluster}-{clade}.gff"
    output:
        "analysis/locate_viruses/segments/{cluster}-{clade}.gff3"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/simplify_gff.py"

checkpoint gff3_to_fna:
    input:
        gff = "analysis/locate_viruses/segments/{cluster}-{clade}.gff3",
        fna = "analysis/locate_viruses/clusters/{cluster}.fna"
    output:
        "analysis/locate_viruses/segments/{cluster}-{clade}.fna"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/gff_fna.py"

rule locate_get_lens:
    input:
        "analysis/locate_viruses/segments/{cluster}-{clade}.fna"
    output:
        "analysis/locate_viruses/segments/{cluster}-{clade}.lens.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit fx2tab -nil {input} | awk '{{print$1,$1,$2}}' > {output}"

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
    params:
        dist = 400
    shell:
        """
        sort -k1,1 -k2,2n {input} | bedtools merge -s -d {params.dist} -c 4,5,6 -o collapse,collapse,distinct | awk '{{split($4,s,",");split($5,c,",");S=0;for(i in s) if(s[i]>S){{S=s[i];C=c[i]}};print$1,$2,$3,S,C,$6}}' OFS=\\\\t > {output}
        """

rule scaf_lens:
    input:
        "analysis/locate_viruses/clusters/cymbo_MaSURCA_assembly_2020.fna"
    output:
        "analysis/locate_viruses/clusters/cymbo_MaSURCA_assembly_2020.lens"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit fx2tab -nl -o {output} {input}"

rule gbk_to_gff3:
    input:
        "../genome_annotation/output/submission/annot_fixed.gbf"
    output:
        "analysis/genome/annot_fixed.gff3"
    conda:
        "envs/bioperl.yaml"
    shell:
        "bp_genbank2gff3.pl -in stdin -out stdout < {input} | perl -pl -MURI::Escape -e '$_=uri_unescape($_)' > {output}"

rule complement:
    input:
        lens = "analysis/locate_viruses/clusters/cymbo_MaSURCA_assembly_2020.lens",
        gff3 = expand("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3", clade = clades)
    output:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed"
    conda:
        "envs/tools.yaml"
    shell:
        "cat {input.gff3} | bedtools sort | bedtools complement -g {input.lens} -i - > {output}"

rule fna_gc:
    input:
        gff3 = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3",
        fna = "analysis/locate_viruses/clusters/cymbo_MaSURCA_assembly_2020.fna"
    output:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3.gc"
    conda:
        "envs/tools.yaml"
    shell:
        "bedtools sort -i {input.gff3} | bedtools merge | bedtools getfasta -fi {input.fna} -bed - | seqkit replace -p '[:-].*' | seqkit fx2tab | bedtools groupby -g 1 -c 2 -o collapse | seqkit tab2fx | seqkit fx2tab -gln -o {output}"

rule complement_gc:
    input:
        bed = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed",
        fna = "analysis/locate_viruses/clusters/cymbo_MaSURCA_assembly_2020.fna"
    output:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed.gc"
    conda:
        "envs/tools.yaml"
    shell:
        "bedtools sort -i {input.bed} | bedtools getfasta -fi {input.fna} -bed - | seqkit replace -p '[:-].*' | seqkit fx2tab | bedtools groupby -g 1 -c 2 -o collapse | seqkit tab2fx | seqkit fx2tab -gln -o {output}"

rule complement_intron_mRNA:
    input:
        annots = "analysis/genome/annot_fixed.gff3",
        compl = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed"
    output:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.mRNA.tab"
    conda:
        "envs/tools.yaml"
    shell:
        "gffread -FCG {input.annots} | bedtools intersect -a - -b {input.compl} | csvgrep --no-header-row -t -c 3 -m mRNA | csvcut -c1,4,5,9 | sed 1d > {output}"

rule count_viruses:
    input:
        MCPs = "analysis/locate_viruses/all_MCPs/{cluster}.bed",
        gff3 = expand("analysis/locate_viruses/segments/{{cluster}}-{clade}.gff3", clade = clades)
    output:
        "output/locate_viruses_{cluster}.svg"
    params:
        clades = clades,
        complete_score = 0.9,
        complete_lens = { "NCLDV": 100000 }
    conda:
        "envs/r.yaml"
    script:
        "scripts/count_viruses.R"

rule scaffold_coverage:
    input:
        "coverage/depth.PLY262_illumina_001_v_PLY262_masurca_assembly.bam.out"
    output:
        "analysis/coverage/{scaffold}.tab"
    shell:
        "grep {wildcards.scaffold} {input} > {output}"

rule plot_coverage:
    input:
        complement_gc = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed.gc",
        clade_gc = expand("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3.gc", clade = clades),
        coverage = "coverage/PLY262_illumina_001_v_PLY262_masurca_assembly.bam.coverage",
        colors = "metadata/clade_colors.tab"
    output:
        "output/total_coverage.pdf",
    params:
        clades = clades,
        scaf_min = 10000
    conda:
        "envs/r.yaml"
    script:
        "scripts/coverage.R"

rule dload_eggnog:
    output:
        "analysis/eggnog/data/eggnog.db"
    conda:
        "envs/eggnog.yaml"
    shell:
        "download_eggnog_data.py -y --data_dir $(dirname {output})"

rule get_scaffold_to_annotate:
    input:
        "clusters/cymbo_MaSURCA_assembly_2020.fasta"
    output:
        "analysis/scaffolds_to_annotate/{scaffold}.fna"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit grep -p {wildcards.scaffold} -o {output} {input}"

rule eggnog:
    input:
        db = "analysis/eggnog/data/eggnog.db",
        faa = "analysis/scaffolds_to_annotate/{genome}.combined.faa"
    output:
        "analysis/scaffolds_to_annotate/{genome}/eggnog.emapper.annotations"
    conda:
        "envs/eggnog.yaml"
    threads:
        20
    shell:
        "emapper.py --cpu {threads} -o eggnog -i {input.faa} --data_dir $(dirname {input.db}) -m diamond --output_dir $(dirname {output})"

rule scaffold_gvog:
    input:
        hmm = "databases/GVDB/output/GVOG.hmm",
        faa = "analysis/scaffolds_to_annotate/{scaffold}.combined.faa"
    output:
        "analysis/scaffolds_to_annotate/{scaffold}_gvog.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmscan -o /dev/null --tblout {output} {input.hmm} {input.faa}"

rule scaffold_interproscan:
    input:
        "analysis/scaffolds_to_annotate/{scaffold}.combined.faa"
    output:
        "analysis/scaffolds_to_annotate/{scaffold}.combined.faa.tsv"
    threads:
        20
    shell:
        "interproscan.sh --disable-precalc -i {input} -d $(dirname {output}) -cpu {threads}"

rule plot_scaffold:
    input:
        eggnog = "analysis/scaffolds_to_annotate/{scaffold}/eggnog.emapper.annotations",
        genes = "analysis/scaffolds_to_annotate/{scaffold}.combined.gff",
        cogs = "metadata/cogs.tab",
        fna = "analysis/scaffolds_to_annotate/{scaffold}.fna",
        cov = "analysis/coverage/{scaffold}.tab",
        cov_total = "coverage/PLY262_illumina_001_v_PLY262_masurca_assembly.bam.coverage",
        viruses = expand("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3", clade = clades),
        colors = "metadata/clade_colors.tab",
        gvog = "analysis/scaffolds_to_annotate/{scaffold}_gvog.txt",
        markers = "metadata/NCLDV_markers.txt",
        iprscan = "analysis/scaffolds_to_annotate/{scaffold}.combined.faa.tsv",
        custom_domains = "metadata/custom_domains.tab"
    output:
        "output/scaffold_{scaffold}.svg"
    params:
        clades = clades
    conda:
        "envs/gmoviz.yaml"
    script:
        "scripts/plot_scaffold.R"

rule summarize_scaffolds:
    input:
        mRNA = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.mRNA.tab",
        complement_gc = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020.complement.bed.gc",
        clade_gc = expand("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3.gc", clade = clades),
        coverage = "coverage/PLY262_illumina_001_v_PLY262_masurca_assembly.bam.coverage",
        MCP = "analysis/locate_viruses/all_MCPs/cymbo_MaSURCA_assembly_2020.bed",
        gff3 = expand("analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3", clade = clades),
        proteinortho = "../genome_annotation/analysis/proteinortho/proteinortho.proteinortho.tsv"
    output:
        "output/viral_scaffolds.tsv"
    params:
        clades = clades,
        min_species = 4
    conda:
        "envs/r.yaml"
    script:
        "scripts/scaffold_summary.R"

rule stringtie_transcripts_intersect:
    input:
        assembly = "../genome_annotation/analysis/stringtie/assembly.gtf",
        virus = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.gff3"
    output:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.stringtie.gff"
    params:
        m = 200
    conda:
        "envs/tools.yaml"
    shell:
        "bedtools intersect -a {input.assembly} -b {input.virus} -f 1 -wb > {output}"

rule stringtie_transcripts_sequence:
    input:
        gff = "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.stringtie.gff",
        genome = "analysis/locate_viruses/clusters/cymbo_MaSURCA_assembly_2020.fna"
    output:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.stringtie.ffn"
    params:
        m = 200
    conda:
        "envs/tools.yaml"
    shell:
        "cut -f1-9 {input.gff} | gffread -Ww- -g {input.genome} | seqkit seq -gm {params.m} > {output}"

rule stringtie_transdecoder:
    input:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.stringtie.ffn"
    output:
        "analysis/locate_viruses/segments/cymbo_MaSURCA_assembly_2020-{clade}.stringtie.faa"
    shadow:
        "minimal"
    conda:
        "envs/transdecoder.yaml"
    shell:
        """
        TransDecoder.LongOrfs -t {input} -O ./ || true
        [ -s longest_orfs.pep ] && mv longest_orfs.pep {output} || touch {output}
        """
