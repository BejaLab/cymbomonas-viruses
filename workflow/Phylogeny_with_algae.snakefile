
rule mcp_algae:
    input:
        expand("output/{CP}.svg", CP = "MCP_PLV"),

rule algae_getorf:
    input:
        "genomes_algae/{genome}.fna"
    output:
        "analysis/getorf_algae/{genome}.faa"
    params:
        minsize = 100,
        maxsize = 100000
    conda:
        "envs/emboss.yaml"
    shell:
        "getorf -minsize {params.minsize} -maxsize {params.maxsize} -filter {input} > {output}"

rule cluster_core_genes:
    input:
        "metadata/queries/{CP}.faa"
    output:
        "analysis/queries/{CP}.cdhit"
    params:
        c = 0.9
    conda:
        "envs/tools.yaml"
    shell:
        "cdhit -i {input} -o {output} -c {params.c} -d 0"

rule CP_align:
    input:
        "metadata/queries/{CP}.faa"
    output:
        "analysis/queries/{CP}.align"
    conda:
        "envs/tools.yaml"
    shell:
        "mafft {input} > {output}"

rule CP_hmmbuild:
    input:
        "analysis/queries/{CP}.align"
    output:
        "analysis/queries/{CP}.hmm"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmbuild {output} {input}"

rule algae_hmmsearch:
    input:
        faa = "analysis/getorf_algae/{genome}.faa",
        hmm = "analysis/queries/{CP}.hmm"
    output:
        "analysis/phylogeny/algae/hmm/{genome}-{CP}.tblout"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmsearch -o /dev/null --tblout {output} {input.hmm} {input.faa}"

rule PLV_hmmsearch:
    input:
        faa = "analysis/PLVs/{genome}.faa",
        hmm = "analysis/queries/{CP}.hmm"
    output:
        "analysis/phylogeny/PLVs/hmm/{genome}-{CP}.tblout"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmsearch -o /dev/null --tblout {output} {input.hmm} {input.faa}"

rule algae_blast:
    input:
        faa = "analysis/getorf_algae/{genome}.faa",
        pdb = "analysis/getorf_algae/{genome}.faa.pdb",
        query = "metadata/queries/{CP}.faa"
    output:
        "analysis/phylogeny/algae/blast/{genome}-{CP}.blast"
    params:
        headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    conda:
        "envs/tools.yaml"
    shell:
        "blastp -query {input.query} -db {input.faa} -outfmt '6 {params.headers}' -out {output}"

rule PLV_blast:
    input:
        faa = "analysis/PLVs/{genome}.faa",
        pdb = "analysis/PLVs/{genome}.faa.pdb",
        query = "metadata/queries/{CP}.faa"
    output:
        "analysis/phylogeny/PLVs/blast/{genome}-{CP}.blast"
    params:
        headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    conda:
        "envs/tools.yaml"
    shell:
        "blastp -query {input.query} -db {input.faa} -outfmt '6 {params.headers}' -out {output}"

rule hmmsearch_extract_algae:
    input:
        fasta = "analysis/getorf_algae/{genome}.faa",
        fai = "analysis/getorf_algae/{genome}.faa.seqkit.fai",
        tblout = "analysis/phylogeny/algae/hmm/{genome}-{CP}.tblout"
    output:
        "analysis/phylogeny/algae/hmm/{genome}-{CP}.faa"
    params:
        evalue = lambda w: evalues[w.CP]
    conda:
        "envs/tools.yaml"
    shell:
        "awk -ve={params.evalue} '!/^#/&&$5<e{{print$1}}' {input.tblout} | sort -u | xargs -r seqkit faidx -f {input.fasta} | seqkit replace -p '^([^ ]+)' -r '$1 {wildcards.genome}' -o {output}"

rule hmmsearch_extract_PLV:
    input:
        fasta = "analysis/PLVs/{genome}.faa",
        fai = "analysis/PLVs/{genome}.faa.seqkit.fai",
        tblout = "analysis/phylogeny/PLVs/hmm/{genome}-{CP}.tblout"
    output:
        "analysis/phylogeny/PLVs/hmm/{genome}-{CP}.faa"
    params:
        evalue = lambda w: evalues[w.CP]
    conda:
        "envs/tools.yaml"
    shell:
        "awk -ve={params.evalue} '!/^#/&&$5<e{{print$1}}' {input.tblout} | sort -u | xargs -r seqkit faidx -f {input.fasta} | seqkit replace -p '^([^ ]+)' -r '$1 {wildcards.genome}' -o {output}"

rule blast_extract_algae:
    input:
        fasta = "analysis/getorf_algae/{genome}.faa",
        fai = "analysis/getorf_algae/{genome}.faa.seqkit.fai",
        blast = "analysis/phylogeny/algae/blast/{genome}-{CP}.blast"
    output:
        "analysis/phylogeny/algae/blast/{genome}-{CP}.faa"
    params:
        evalue = 1e-10
    conda:
        "envs/tools.yaml"
    shell:
        "awk -ve={params.evalue} '$11<e' {input.blast} | cut -f2 | sort -u | xargs -r seqkit faidx -f {input.fasta} > {output}"

rule blast_extract_PLVs:
    input:
        fasta = "analysis/PLVs/{genome}.faa",
        fai = "analysis/PLVs/{genome}.faa.seqkit.fai",
        blast = "analysis/phylogeny/PLVs/blast/{genome}-{CP}.blast"
    output:
        "analysis/phylogeny/PLVs/blast/{genome}-{CP}.faa"
    params:
        evalue = 1e-10
    conda:
        "envs/tools.yaml"
    shell:
        "awk -ve={params.evalue} '$11<e' {input.blast} | cut -f2 | sort -u | xargs -r seqkit faidx -f {input.fasta} > {output}"

rule algae_hmm:
    input:
        fasta = "analysis/getorf_algae/{genome}.faa",
        hmm = "hmm_algae/{hmm}.hmm"
    output:
        "analysis/hmm_algae/{genome}-{hmm}.tblout"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmsearch -o /dev/null --tblout {output} {input.hmm} {input.fasta}"

#rule blast_extract_cat:
#    input:
#        expand("analysis/phylogeny/algae/blast/{genome}-MCP_NCLDV.faa", genome = algal_genomes),
#        expand("analysis/phylogeny/PLVs/blast/{genome}-MCP_NCLDV.faa", genome = all_genomes)
#    output:
#        "analysis/phylogeny/MCP_NCLDV.fasta"
#    shell:
#        "seqkit seq -gm200 {input} -o {output}"

rule hmm_extract_cat_PLV:
    input:
        "metadata/queries/MCP_PLV_outgroups.faa",
        "metadata/MCP_PLV_metatranscriptomes.faa",
        expand("analysis/phylogeny/algae/hmm/{genome}-MCP_PLV.faa", genome = algal_genomes),
        expand("analysis/phylogeny/PLVs/hmm/{genome}-MCP_PLV.faa", genome = virus_names)
    output:
        "analysis/phylogeny/MCP_PLV.fasta"
    params:
        m = 200
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -gm{params.m} {input} | seqkit rmdup -o {output}"

rule hmm_extract_cat_NCLDV:
    input:
        "metadata/queries/MCP_NCLDV_outgroups.faa",
        expand("analysis/phylogeny/algae/hmm/{genome}-MCP_NCLDV.faa", genome = [ "Phaglo", "Pharex", "Phaant" ]),
        expand("analysis/phylogeny/PLVs/hmm/{genome}-MCP_NCLDV.faa", genome = virus_names)
    output:
        "analysis/phylogeny/MCP_NCLDV.fasta"
    params:
        m = 200
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -gm{params.m} {input} | seqkit rmdup -o {output}"

rule fasta_hmmalign:
    input:
        fasta = "analysis/phylogeny/{profile}.fasta",
        hmm = "metadata/queries/{profile}_hmmalign.hmm"
    output:
        "analysis/phylogeny/{profile}.a2m"
    conda:
        "envs/tools.yaml"
    shell:
        "hmmalign --trim --outformat A2M {input.hmm} {input.fasta} > {output}"

rule hmmalign_trim:
    input:
        a2m = "analysis/phylogeny/{profile}.a2m",
        fai = "analysis/phylogeny/{profile}.a2m.seqkit.fai"
    output:
        "analysis/phylogeny/{profile}.a2m.trim"
    params:
        m = 300
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit replace -sp [a-z] {input.a2m} | seqkit seq -nigm{params.m} | xargs seqkit faidx {input.a2m} | seqkit replace -sp [a-z] | seqkit rmdup -so {output}"

rule hmmalign_select_short:
    input:
        a2m = "analysis/phylogeny/{profile}.a2m",
        trim = "analysis/phylogeny/{profile}.a2m.trim"
    output:
        "analysis/phylogeny/{profile}.a2m.trim.short"
    params:
        M = 299
    shell:
        "seqkit replace -sp [a-z] {input.a2m} | seqkit seq -nigM{params.M} | xargs seqkit faidx {input.a2m} | seqkit replace -sp [a-z] | seqkit rmdup -so {output}"

rule raxml_evaluate:
    input:
        fasta  = "analysis/phylogeny/{profile}.a2m.trim",
        iqtree = "analysis/phylogeny/{profile}.treefile"
    output:
        "analysis/phylogeny/{profile}.raxml.bestTree",
        "analysis/phylogeny/{profile}.raxml.bestModel"
    params:
        prefix = "analysis/phylogeny/{profile}",
        seed = 123,
        model = lambda w: models[w.profile]
    log:
        "analysis/phylogeny/{profile}.raxml.log"
    conda:
        "envs/placement.yaml"
    shell:
        "raxml-ng --redo --evaluate --msa {input.fasta} --tree {input.iqtree} --model {params.model} --seed {params.seed} --prefix {params.prefix} &> {log}"

rule epa:
    input:
        fasta = "analysis/phylogeny/{profile}.a2m.trim",
        tree  = "analysis/phylogeny/{profile}.raxml.bestTree",
        model = "analysis/phylogeny/{profile}.raxml.bestModel",
        short = "analysis/phylogeny/{profile}.a2m.trim.short"
    output:
        "analysis/phylogeny/{profile}_epa/epa_result.jplace"
    params:
        dir = "analysis/phylogeny/{profile}_epa"
    log:
        "analysis/phylogeny/{profile}_epa/epa.log"
    conda:
        "envs/placement.yaml"
    shell:
        "epa-ng --redo -s {input.fasta} -t {input.tree} --model {input.model} -q {input.short} -w {params.dir} &> {log}"

rule gappa:
    input:
        "analysis/phylogeny/{profile}_epa/epa_result.jplace"
    output:
        "analysis/phylogeny/{profile}_epa/epa_result.newick"
    params:
        out_dir = "analysis/phylogeny/{profile}_epa"
    conda:
        "envs/placement.yaml"
    shell:
        "gappa examine graft --jplace-path {input} --fully-resolve --out-dir {params.out_dir}"

rule algae_iqtree_hmm:
    input:
        "analysis/phylogeny/{profile}.a2m.trim"
        # "{prefix}.trimal.uniq"
    output:
        "analysis/phylogeny/{profile}.treefile"
    params:
        seed = 123,
        B = 1000,
        prefix = "analysis/phylogeny/{profile}"
    threads:
        4
    conda:
        "envs/tools.yaml"
    shell:
        "iqtree2 -s {input} --prefix {params.prefix} -redo --alrt 1000 -B 1000 --seed {params.seed} -T {threads}"

rule algae_ggtree_MCP_NCLDV:
    input:
        tree = "analysis/phylogeny/MCP_NCLDV.treefile",
        # tree = "analysis/phylogeny/MCP_NCLDV.treefile",
        fasta = "analysis/phylogeny/MCP_NCLDV.fasta",
        synonyms = "metadata/organisms.txt",
        outgroups = "metadata/queries/MCP_NCLDV_outgroups.faa"
    output:
        image = "output/MCP_NCLDV.svg",
        jtree = "output/MCP_NCLDV.jtree"
    params:
        outgroup_rooting = False
    conda:
        "envs/r.yaml"
    script:
        "scripts/ggtree.R"

rule algae_ggtree_MCP_PLV:
    input:
        tree = "analysis/phylogeny/MCP_PLV_epa/epa_result.newick",
        fasta = "analysis/phylogeny/MCP_PLV.fasta",
        synonyms = "metadata/organisms.txt",
        outgroups = "metadata/queries/MCP_PLV_outgroups.faa"
    output:
        image = "output/MCP_PLV.svg",
        jtree = "output/MCP_PLV.jtree"
    params:
        outgroup_rooting = True
    conda:
        "envs/r.yaml"
    script:
        "scripts/ggtree.R"

rule algae_ggtree_MCP:
    input:
        tree = "analysis/phylogeny/MCP.treefile",
        fasta = expand("analysis/blast_algae/{genome}-MCP.faa", genome = algal_genomes),
        tblout = expand("analysis/hmm_algae/{genome}-{hmm}.tblout", genome = algal_genomes, hmm = hmm_algae),
        synonyms = "metadata/organisms.txt",
        hmm = expand("hmm_algae/{hmm}.hmm", hmm = hmm_algae)
    output:
        "analysis/phylogeny/MCP.svg"
    conda:
        "envs/r.yaml"
    script:
        "scripts/ggtree.R"
