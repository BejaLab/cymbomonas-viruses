
rule viruses_gvmags_faa:
    input:
        "databases/GVMAGs/data"
    output:
        "analysis/PLV_annotate/gvmags-{virus}.faa"
    shell:
        "cat {input}/GVMAGs_*/GVMAGs_*_faa/{wildcards.virus}.faa | tr '|' '_' > {output}"

rule viruses_moniruzzaman_faa:
    input:
        "databases/Moniruzzaman20/prot_for_posting/{virus}.fa.dc.faa"
    output:
        "analysis/PLV_annotate/moniruzzaman20-{virus}.faa"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit replace -p ^ -r {wildcards.virus}_ {input} | seqkit seq -i | seqkit replace -sp '[*]' -o {output}"

rule viruses_repbase_fna:
    input:
        "databases/RMRB/Libraries/RMRBSeqs.embl.fasta"
    output:
        "analysis/PLV_annotate/prokka-{genome}.fna"
    wildcard_constraints:
        genome = '|'.join([ v for v in viruses if viruses[v]['dataset'] == 'repbase' ])
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit grep -p {wildcards.genome} -o {output} {input}"

rule viruses_ncbi_fna:
    output:
        "analysis/PLV_annotate/prokka-{genome}.fna"
    wildcard_constraints:
        genome = '|'.join([ v for v in viruses if viruses[v]['dataset'] == 'ncbi' ])
    params:
        id = lambda w: viruses[w.genome]['accession']
    resources:
        ncbi = 1
    conda:
        "envs/tools.yaml"
    shell:
        "efetch -db nuccore -id {params.id} -format fasta > {output}"

rule viruses_ncbi_gbk:
    output:
        "analysis/PLV_annotate/ncbi-{genome}.gbk"
    params:
        id = lambda w: viruses[w.genome]['accession']
    resources:
        ncbi = 1
    conda:
        "envs/tools.yaml"
    shell:
        "efetch -db nuccore -id {params.id} -format gb > {output}"

rule cp_gbk:
    input:
        lambda w: expand("analysis/PLV_annotate/{annotations}-{genome}.gbk", annotations = viruses[w.genome]['annotations'], genome = w.genome)
    output:
        "analysis/PLVs/{genome}.gbk"
    shell:
        "cp {input} {output}"

rule cp_lens:
    input:
        lambda w: expand("analysis/PLV_annotate/{annotations}-{genome}.lens.txt", annotations = viruses[w.genome]['annotations'], genome = w.genome)
    output:
        "analysis/PLVs/{genome}.lens.txt"
    shell:
        "cp {input} {output}"

rule mv_faa:
    input:
        lambda w: expand("analysis/PLV_annotate/{annotations}-{genome}.faa", annotations = viruses[w.genome]['annotations'], genome = w.genome)
    output:
        "analysis/PLVs/{genome}.faa"
    shell:
        "mv {input} {output}"

rule prokka_gbk:
    input:
        "analysis/PLV_annotate/prokka-{genome}.fna"
    output:
        "analysis/PLV_annotate/prokka-{genome}.gbk"
    params:
        outdir = "analysis/PLVs/ncbi_annotate/",
        locustag = lambda w: re.sub('[^A-Z]', '', w.genome.upper())
    threads:
        8
    # NB: prokka cannot be installed with conda
    # conda:
    #    "envs/prokka.yaml"
    shell:
        """
        prokka --force --cpus {threads} --prefix {wildcards.genome} --locustag {params.locustag} --kingdom Viruses --outdir {params.outdir}/{wildcards.genome} {input}
        cp {params.outdir}/{wildcards.genome}/{wildcards.genome}.gbk {output}
        """

rule extract_faa_bellas:
    input:
        "databases/Bellas_Sommaruga/input/All_proteins.faa"
    output:
        "analysis/PLV_annotate/bellas-{virus}.faa"
    params:
        search = lambda w: viruses[w.virus]['search']
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit grep -rp {params.search} {input} | seqkit replace -p '$' -r ' {wildcards.virus}' -o {output}"

rule extract_gbk_manual:
    input:
        "databases/manual_viruses/{virus}.gbk"
    output:
        "analysis/PLV_annotate/manual-{virus}.gbk"
    shell:
        "cp {input} {output}"

rule extract_faa_manual:
    input:
        "analysis/PLV_annotate/manual-{virus}.gbk"
    output:
        "analysis/PLV_annotate/manual-{virus}.faa"
    conda:
        "envs/bioperl.yaml"
    shell:
        "perl workflow/helpers/biotags.pl -i {input} -p CDS -T id -t locus_tag,spliced_seq | awk '{{printf\"%s %s\\t%s\\n\",$2,$1,$3}}' | seqkit tab2fx | seqkit translate --trim -o {output}"

rule extract_gbk_yutin:
    input:
        "databases/Yutin2015.gbk"
    output:
        "analysis/PLV_annotate/yutin-{virus}.gbk"
    conda:
        "envs/tools.yaml"
    shell:
        "awk -vl={wildcards.virus} '/^LOCUS/{{p=$2==l}}p' {input} | sed -E 's/\\x93|\\x94/\"/g' > {output}"

rule extract_faa_yutin:
    input:
        "analysis/PLV_annotate/yutin-{virus}.gbk"
    output:
        "analysis/PLV_annotate/yutin-{virus}.faa"
    conda:
        "envs/bioperl.yaml"
    shell:
        "perl workflow/helpers/biotags.pl -i {input} -p CDS -T id -t locus_tag,translation | awk '{{printf\"%s %s\\t%s\\n\",$2,$1,$3}}' | seqkit tab2fx -o {output}"

rule extract_faa_ncbi_annotate:
    input:
        "analysis/PLV_annotate/prokka-{virus}.gbk"
    output:
        "analysis/PLV_annotate/prokka-{virus}.faa"
    conda:
        "envs/bioperl.yaml"
    shell:
        "perl workflow/helpers/biotags.pl -i {input} -p CDS -T description -t locus_tag,translation | awk -F\\\\t '$2{{printf\"%s %s\\t%s\\n\",$2,$1,$3}}' | seqkit tab2fx | seqkit seq -go {output}"

rule extract_faa_ncbi:
    input:
        "analysis/PLV_annotate/ncbi-{virus}.gbk"
    output:
        "analysis/PLV_annotate/ncbi-{virus}.faa"
    conda:
        "envs/bioperl.yaml"
    shell:
        "perl workflow/helpers/biotags.pl -i {input} -p CDS -T description -t protein_id,translation | awk -F\\\\t '$2{{printf\"%s %s\\t%s\\n\",$2,$1,$3}}' | seqkit tab2fx -o {output}"

rule get_lens_gbk:
    input:
        "analysis/PLV_annotate/{source}-{virus}.gbk",
    wildcard_constraints:
        source = "ncbi|prokka|yutin|manual"
    output:
        "analysis/PLV_annotate/{source}-{virus}.lens.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "grep ^LOCUS {input} | awk -vV={wildcards.virus} '{{print V,$2,$3}}' > {output}"

rule get_lens_moniruzzaman:
    input:
        "databases/Moniruzzaman20/final_bins_nucl/{virus}.dc.fa"
    output:
        "analysis/PLV_annotate/moniruzzaman20-{virus}.lens.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit fx2tab -nil {input} | awk -vV={wildcards.virus} '{{print V,$1,$2}}' > {output}"

rule get_lens_bellas:
    input:
        "databases/Bellas_Sommaruga/input/All_proteins.faa"
    output:
        "analysis/PLV_annotate/bellas-{virus}.faa"
    params:
        search = lambda w: viruses[w.virus]['search']
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit grep -rp {params.search} {input} | seqkit seq -ni | awk -F_[-]_ '{{c=$1}}$2~/_len/{{c=$2}}{{print c}}' | sort -u | awk -F_ -vV={wildcards.virus} '{{for(i=1;i<NF;i++)if($i~/^len/){{print V,$0,$(i+1);next}}}} > {output}"

rule get_lens_gvmags:
    input:
        "databases/GVMAGs/data"
    output:
        "analysis/PLV_annotate/gvmags-{virus}.lens.txt"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit fx2tab -nil {input}/GVMAGs_*/GVMAGs_*_fna/{wildcards.virus}.fna | awk -vV={wildcards.virus} '{{print V,$1,$2}}' > {output}"
