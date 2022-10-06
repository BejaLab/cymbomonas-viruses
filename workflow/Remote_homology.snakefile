
rule remote_homology:
    input:
        expand("analysis/PLV_segments_hh/{segment}-hhsearch-pfam.tsv", segment = segment_names),
        expand("analysis/PLVs_hh/{segment}-hhsearch-pfam.tsv", segment = small_viruses + Mesomimi_viruses),
        expand("output/core_genes_{coverage}.pdf", coverage = [ 40, 60 ]),

rule collect_proteins:
    input:
        expand("analysis/PLVs/{genome}.faa", genome = small_viruses + Mesomimi_viruses),
        expand("analysis/PLV_segments/{segment}.faa", segment = segment_names)
    output:
        "analysis/hhblits_db/All_proteins.faa"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit seq -o {output} {input}"

rule segments_faa_ffindex:
    input:
        "analysis/PLVs/{segment}.faa"
    output:
        data  = "analysis/PLVs/{segment}.ffdata",
        index = "analysis/PLVs/{segment}.ffindex"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "ffdb fasta -d {output.data} -i {output.index} {input}"

rule PLV_faa_ffindex:
    input:
        "analysis/PLV_segments/{segment}.faa"
    output:
        data  = "analysis/PLV_segments/{segment}.ffdata",
        index = "analysis/PLV_segments/{segment}.ffindex"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "ffdb fasta -d {output.data} -i {output.index} {input}"

rule ffindex_num:
    input:
        "{prefix}.ffindex"
    output:
        "{prefix}.ffindex.num"
    shell:
        "cut -f2,3 {input} | nl -w1 > {output}"

rule hhblits_a3m:
    input:
        data  = "analysis/{dir}/{segment}.ffdata",
        index = "analysis/{dir}/{segment}.ffindex.num",
        msa_ffdata = "analysis/hhblits_db/All_proteins_msa_sequence.ffdata",
        msa_cs219_ffdata = "analysis/hhblits_db/All_proteins_msa_cs219.ffdata"
    output:
        data  = "analysis/{dir}_hh/{segment}-hhblits.a3m.ffdata",
        index = "analysis/{dir}_hh/{segment}-hhblits.a3m.ffindex"
    params:
        db = "analysis/hhblits_db/All_proteins_msa",
        n = 3,
        e = 1e-5
    conda:
        "envs/soedinglab.yaml"
    shell:
        "ffindex_apply -q -d {output.data} -i {output.index} {input.data} {input.index} -- hhblits -cpu 1 -v 0 -b 1 -z 1 -d {params.db} -i stdin -oa3m stdout -e {params.e} -n {params.n}"

rule hhblits_ss:
    input:
        data  = "analysis/{dir}_hh/{segment}-hhblits.a3m.ffdata",
        index = "analysis/{dir}_hh/{segment}-hhblits.a3m.ffindex"
    output:
        data  = "analysis/{dir}_hh/{segment}-hhblits-ss.a3m.ffdata",
        index = "analysis/{dir}_hh/{segment}-hhblits-ss.a3m.ffindex"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "ffindex_apply -q -d {output.data} -i {output.index} {input.data} {input.index} -- $CONDA_PREFIX/scripts/addss.pl -v 0 /dev/stdin /dev/stdout"

rule hhsearch_self:
    input:
        data  = "analysis/{dir}_hh/{segment}-hhblits.a3m.ffdata",
        index = "analysis/{dir}_hh/{segment}-hhblits.a3m.ffindex",
        contxt = "databases/hhsuite/data/context_data.crf",
        msa_ffdata = "analysis/hhblits_db/All_proteins_msa_sequence.ffdata",
        msa_cs219_ffdata = "analysis/hhblits_db/All_proteins_msa_cs219.ffdata"
    output:
        data  = "analysis/{dir}_hh/{segment}-hhsearch-self.hhr.ffdata",
        index = "analysis/{dir}_hh/{segment}-hhsearch-self.hhr.ffindex"
    params:
        db = "analysis/hhblits_db/All_proteins_msa",
        BZ  = 250,
        p   = 20,
        ssm = 2
    conda:
        "envs/soedinglab.yaml"
    shell:
        "ffindex_apply -q -d {output.data} -i {output.index} {input.data} {input.index} -- hhsearch -v 1 -cpu 1 -b 1 -z 1 -i stdin -d {params.db} -o stdout -p {params.p} -Z {params.BZ} -B {params.BZ} -ssm {params.ssm} -contxt {input.contxt}"

rule hhsearch_pfam:
    input:
        data  = "analysis/{dir}_hh/{segment}-hhblits-ss.a3m.ffdata",
        index = "analysis/{dir}_hh/{segment}-hhblits-ss.a3m.ffindex",
        contxt = "databases/hhsuite/data/context_data.crf",
        db = "databases/hhsuite/data/pfam_hhm.ffdata"
    output:
        data  = "analysis/{dir}_hh/{segment}-hhsearch-pfam.hhr.ffdata",
        index = "analysis/{dir}_hh/{segment}-hhsearch-pfam.hhr.ffindex"
    params:
        db = "databases/hhsuite/data/pfam",
        BZ  = 250,
        p   = 20,
        ssm = 2
    conda:
        "envs/soedinglab.yaml"
    shell:
        "ffindex_apply -q -d {output.data} -i {output.index} {input.data} {input.index} -- hhsearch -v 1 -cpu 1 -b 1 -z 1 -i stdin -d {params.db} -o stdout -p {params.p} -Z {params.BZ} -B {params.BZ} -ssm {params.ssm} -contxt {input.contxt}"

rule hhsuite_parse:
    input:
        data  = "{prefix}.hhr.ffdata",
        index = "{prefix}.hhr.ffindex"
    output:
        "{prefix}.tsv"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "ffindex_apply -q {input.data} {input.index} -- Rscript workflow/helpers/parse_hhsuite.R | sed '2,${{/^Query/d}}' > {output}"

rule protein_createdb_bellas:
    input:
        "databases/Bellas_Sommaruga/input/All_proteins.faa"
    output:
        "analysis/hhblits_db/Bellas"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "mmseqs createdb {input} {output}"

rule protein_createdb:
    input:
        "analysis/hhblits_db/All_proteins.faa"
    output:
        "analysis/hhblits_db/All_proteins"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "mmseqs createdb {input} {output}"

rule protein_cluster:
    input:
        "analysis/hhblits_db/All_proteins"
    output:
        "analysis/hhblits_db/All_proteins_clu.0",
        "analysis/hhblits_db/All_proteins_clu.dbtype",
        "analysis/hhblits_db/All_proteins_clu.index"
    params:
        output_prefix = "analysis/hhblits_db/All_proteins_clu",
        tmp = "analysis/hhblits_db/tmp",
        min_seq_id = 0.3,
        coverage = 0.8
    threads:
        10
    conda:
        "envs/soedinglab.yaml"
    shell:
        "mmseqs cluster --min-seq-id {params.min_seq_id} -c {params.coverage} --threads {threads} {input} {params.output_prefix} {params.tmp}"

rule protein_result2msa:
    input:
        clu = "analysis/hhblits_db/All_proteins_clu.0",
        db = "analysis/hhblits_db/All_proteins"
    output:
        expand("analysis/hhblits_db/All_proteins_msa_{type}.{ext}", type = [ "ca3m", "header", "sequence" ], ext = [ "ffdata", "ffindex" ])
    params:
        input_prefix  = "analysis/hhblits_db/All_proteins_clu",
        output_prefix = "analysis/hhblits_db/All_proteins_msa"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "mmseqs result2msa {input.db} {input.db} {params.input_prefix} {params.output_prefix} --msa-format-mode 1"

rule protein_cstranslate:
    input:
        msa_header_ffdata = "analysis/hhblits_db/All_proteins_msa_header.ffdata",
        db = "analysis/hhblits_db/All_proteins",
        cs219_lib = "databases/hhsuite/data/cs219.lib",
        context = "databases/hhsuite/data/context_data.lib"
    output:
        "analysis/hhblits_db/All_proteins_msa_cs219.ffdata",
        "analysis/hhblits_db/All_proteins_msa_cs219.ffindex"
    params:
        input_prefix  = "analysis/hhblits_db/All_proteins_msa",
        output_prefix = "analysis/hhblits_db/All_proteins_msa_cs219"
    threads:
        30
    # NB: conda doesn't have mpi version of hhsuite
    #conda:
    #    "envs/soedinglab.yaml"
    shell:
        "mpirun -np {threads} cstranslate_mpi -i {params.input_prefix} -o {params.output_prefix} -A {input.cs219_lib} -D {input.context} -x 0.3 -c 4 -I ca3m -b"

rule createtsv:
    input:
        db = "analysis/hhblits_db/All_proteins",
        clu = "analysis/hhblits_db/All_proteins_clu.0"
    output:
        "analysis/hhblits_db/All_proteins_clu.tsv"
    params:
        clu_prefix = "analysis/hhblits_db/All_proteins_clu"
    conda:
        "envs/soedinglab.yaml"
    shell:
        "mmseqs createtsv {input.db} {input.db} {params.clu_prefix} {output}"

rule mcl_abc:
    input:
        clu_tsv = "analysis/hhblits_db/All_proteins_clu.tsv",
        segment_tsv = expand("analysis/PLV_segments_hh/{segment}-hhsearch-self.tsv", segment = segment_names),
        virus_tsv = expand("analysis/PLVs_hh/{segment}-hhsearch-self.tsv", segment = small_viruses + Mesomimi_viruses),
    output:
        "analysis/hh_mcl/abc_{coverage}.tsv"
    params:
        coverage = "{coverage}",
        probab = 90,
        coverage_q_frag = 80,
        coverage_t_frag = 25,
        identities_frag = 30,
        probab_frag = 99
    conda:
        "envs/r.yaml"
    script:
        "scripts/abc_graph.R"

rule mcl_run:
    input:
        "analysis/hh_mcl/abc_{coverage}.tsv"
    output:
        "analysis/hh_mcl/mcl_{coverage}.txt"
    params:
        I = 2,
        scheme = 7
    conda:
        "envs/tools.yaml"
    threads:
        10
    shell:
        "mcl {input} -o {output} --abc -scheme {params.scheme} -te {threads} -I {params.I}"

rule mcl_analyze:
    input:
        mcl = "analysis/hh_mcl/mcl_{coverage}.txt",
        segment_ffindex  = expand("analysis/PLV_segments/{segment}.ffindex", segment = segment_names),
        virus_ffindex    = expand("analysis/PLVs/{virus}.ffindex", virus = small_viruses + Mesomimi_viruses),
        segment_pfam     = expand("analysis/PLV_segments_hh/{segment}-hhsearch-pfam.tsv", segment = segment_names),
        virus_pfam       = expand("analysis/PLVs_hh/{virus}-hhsearch-pfam.tsv", virus = small_viruses + Mesomimi_viruses),
        virus_metadata   = "metadata/viruses.tsv",
        segment_metadata = "metadata/viral_segments.tsv",
        genes_cluster    = "metadata/genes_cluster.tsv",
        genes_pfam       = "metadata/genes_pfam.tsv",
        jtree            = "output/MCP_PLV.jtree",
        colors           = "metadata/subclade_colors.txt",
        plv_order        = "metadata/plv_order.txt",
        family_colors    = "metadata/family_colors.txt"
    output:
        core_genes = "output/core_genes_{coverage}.pdf",
        bipartite  = "output/bipartite_{coverage}.pdf",
        data       = "output/mcl_genes_{coverage}.tsv",
        reduced_tree = "output/MCP_PLV_reduced_{coverage}.svg"
    params:
        probab = 80,
        chosen_clades = [ "Gezel", "Dwarf", "Mesomimi" ],
        clade_levels  = [ "Mesomimi", "Dwarf", "Lavidaviridae", "TVS", "Gezel" ],
        ref_genomes   = [ "Sputnik", "Mavirus_Spezl", "TVV_S1", "Dialut-1a", "Dialut-1b", "Dialut-2" ]
    conda:
        "envs/r.yaml"
    script:
        "scripts/mcl_graph.R"

