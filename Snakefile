
genomes, = glob_wildcards("genomes/{genome}.fna")
hmms,    = glob_wildcards("hmm/{hmm}.hmm")

rule all:
	input:
		expand("analysis/hmmsearch/{genome}-{hmm}.txt", genome = genomes, hmm = hmms),
		expand("analysis/segments/{genome}/segments.pfam.txt", genome = genomes),
		expand("analysis/segments/{genome}/segments.gvog.tblout", genome = genomes),
		expand("analysis/segments/{genome}/segments.vogdb.tblout", genome = genomes),
		expand("analysis/segments/{genome}/segments.bellas.tblout", genome = genomes),
		expand("analysis/segments/{genome}/segments.gbk", genome = genomes)

rule cat_hmm:
	input:
		expand("hmm/{hmm}.hmm", hmm = hmms)
	output:
		"analysis/profiles.hmm"
	shell:
		"cat {input} > {output}"

rule getorf:
	input:
		"genomes/{genome}.fna"
	output:
		"analysis/getorf/{genome}.faa",
	shell:
		"getorf -minsize 100 -filter {input} > {output}"

rule hmmsearch:
	input:
		fasta = "analysis/getorf/{genome}.faa", hmm = "hmm/{hmm}.hmm"
	output:
		"analysis/hmmsearch/{genome}-{hmm}.txt"
	shell:
		"hmmsearch --max -o /dev/null --tblout {output} {input.hmm} {input.fasta}"

rule parse_hmmsearch:
	input:
		expand("analysis/hmmsearch/{{genome}}-{hmm}.txt", hmm = hmms)
	output:
		figure = "images/segments_{genome}.pdf",
		bed = "analysis/segments/{genome}/segments.bed"
	params:
		e_value_threshold = 0.001,
		distance_threshold = 100000,
		genes_threshold = 2,
	script:
		"scripts/parse_hmmsearch.R"

rule segments_extract:
	input:
		genome = "genomes/{genome}.fna", bed = "analysis/segments/{genome}/segments.bed"
	output:
		"analysis/segments/{genome}/segments.fna"
	params:
		flanks = 50000
	shell:
		"seqkit subseq --up-stream {params.flanks} --down-stream {params.flanks} --bed {input.bed} {input.genome} | cut -f1 -d: > {output}"

rule segments_faidx:
	input:
		"analysis/segments/{genome}/segments.fna"
	output:
		"analysis/segments/{genome}/segments.fna.fai"
	shell:
		"seqkit faidx {input}"

rule segments_genemarks:
	input:
		"analysis/segments/{genome}/segments.fna"
	output:
		"analysis/segments/{genome}/segments.fna.gff3"
	params:
		dirname = "analysis/segments"
	shell:
		"""
		cd {params.dirname}/{wildcards.genome}
		gmsn.pl --format GFF3 --virus segments.fna
		"""

rule segments_gffread:
	input:
		fna = "analysis/segments/{genome}/segments.fna",
		gff3 = "analysis/segments/{genome}/segments.fna.gff3"
	output:
		"analysis/segments/{genome}/segments.faa"
	shell:
		"gffread -g {input.fna} -y {output} {input.gff3}"

rule segments_gvog:
	input:
		hmm = "databases/GVDB/output/GVOG.hmm",
		fasta = "analysis/segments/{genome}/segments.faa"
	output:
		tblout = "analysis/segments/{genome}/segments.gvog.tblout",
		domtblout = "analysis/segments/{genome}/segments.gvog.domtblout"
	threads:
		4
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
	shell:
		"hmmscan --cpu {threads} -o /dev/null --tblout {output.tblout} --domtblout {output.domtblout} {input.hmm} {input.fasta}"

rule segments_pfam:
	input:
		"analysis/segments/{genome}/segments.faa"
	output:
		"analysis/segments/{genome}/segments.pfam.txt"
	params:
		pfam_dir = "databases/Pfam"
	threads:
		4
	shell:
		"pfam_scan.pl -fasta {input} -dir {params.pfam_dir} -outfile {output} -cpu {threads}"

rule segments_genbank:
	input:
		fna = "analysis/segments/{genome}/segments.fna",
		gff3 = "analysis/segments/{genome}/segments.fna.gff3",
		pfam_db = "databases/Pfam/Pfam-A.hmm",
		tblout = [
			"analysis/segments/{genome}/segments.gvog.tblout",
			"analysis/segments/{genome}/segments.vogdb.tblout",
			"analysis/segments/{genome}/segments.bellas.tblout"
		],
		pfam = "analysis/segments/{genome}/segments.pfam.txt"
	output:
		"analysis/segments/{genome}/segments.gbk"
	params:
		evalue = 1e-5
	script:
		"scripts/gff2gbk.py"

rule segments_split:
	input:
		"analysis/segments/{genome}/segments.faa"
	output:
		directory("analysis/segments/{genome}/split")
	shell:
		"seqkit split -i -O {output} {input}"

rule segments_hhblits:
	input:
		"analysis/segments/{genome}/split"
	output:
		directory("analysis/segments/{genome}/hhblits")
	params:
		db = "databases/hhsuite/UniRef30_2020_06",
		n = 2,
		e = 1e-3,
		BZ = 250
	threads:
		workflow.cores
	shell:
		"""
		find {input} -name '*.faa' | \
			parallel -j{threads} hhblits -cpu 1 -v 0 -b 1 -z 1 -i {{}} -d {params.db} -o /dev/null -oa3m stdout -e {params.e} -n {params.n} -B {params.BZ} -Z {params.BZ} \| \
			addss.pl -i stdin -o {output}/{{/.}}.a3m
		"""

rule segments_hhsearch:
	input:
		"analysis/segments/{genome}/hhblits"
	output:
		directory("analysis/segments/{genome}/hhsearch")
	params:
		contxt = "databases/hhsuite/context_data.crf",
		db  = "vogdb/vog",
		BZ  = 250,
		p   = 20,
		ssm = 2
	threads:
		workflow.cores
	shell:
		"""
		mkdir -p {output}
		find {input} -name '*.a3m' | \
			parallel -j{threads} hhsearch -cpu 1 -b 1 -z 1 -i {{}} -d {params.db} -o {output}/{{/.}}.hhr -p {params.p} -Z {params.BZ} -B {params.BZ} -ssm {params.ssm} -contxt {params.contxt}
		"""

rule scaffolds_prokka:
	input:
		fasta = "analysis/viral_segments/{genome}.fna", hmm = "analysis/profiles.hmm"
	output:
		"analysis/prokka/{genome}.gff3"
	params:
		outdir = "analysis/prokka"
	shell:
		"prokka --force --kingdom Viruses --hmms {input.hmm} --outdir {params.outdir} --prefix {wildcards.genome} {input.fasta}"
