
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict

ref_annots = []

gff3_file = str(snakemake.input["gff3"])
fna_file  = str(snakemake.input["fna"])

pfam_file = str(snakemake.input["pfam"])
pfam_db = str(snakemake.input["pfam_db"])
tblout_files = snakemake.input["tblout"]

output_file = str(snakemake.output)

evalue_threshold = snakemake.params["evalue"]

aliases = {
	"MCP": "Major capsid protein (MCP)",
	"mCP": "Minor capsid protein (mCP)",
	"uncharacterized protein": "-"
}

descs = {}
with open(pfam_db) as fh:
	for line in fh:
		if line.startswith('ACC'):
			key, ACC = line.rstrip().split()
		if line.startswith('DESC'):
			key, DESC = line.rstrip().split(maxsplit = 1)
			descs[ACC] = DESC

hits = defaultdict(list)

with open(pfam_file) as fh:
	for line in fh:
		if not line.startswith('#') and line != '\n':
			seq_id, aln_start, aln_end, env_start, env_end, hmm_acc, hmm_name, type, hmm_start, hmm_end, hmm_len, bit_score, E_value, significance, clan = line.rstrip().split()
			hit = "%s: %s (%s)" % (hmm_acc, descs[hmm_acc], E_value)
			hits[seq_id].append(hit)

for fname in tblout_files:
	with open(fname) as fh:
		for line in fh:
			if not line.startswith('#'):
				t_name, t_acc, q_name, q_acc, E_value, score, bias, best_E_value, best_score, best_bias, exp, reg, clu, ov, env, dom, rep, inc, desc = line.rstrip().split(maxsplit = 18)
				if desc in aliases:
					desc = aliases[desc]
				definition = ''
				if t_name != desc and desc != '-':
					definition = ' ' + desc
				if float(E_value) <= evalue_threshold:
					hit = "%s:%s (%s)" % (t_name, definition, E_value)
					hits[q_name].append(hit)

features = defaultdict(list)
with open(gff3_file) as fh:
	for line in fh:
		if not line.startswith('#') and line != '\n':
			seqname, source, feature_type, start, end, score, strand, frame, attributes = line.rstrip().split('\t')
			if feature_type == 'CDS':
				atts = {}
				for att in attributes.split(';'):
					key, value = att.split('=')
					atts[key] = value
				parent = atts['Parent']
				loc = FeatureLocation(int(start) - 1, int(end))
				feature = SeqFeature(loc, type = feature_type, strand = int(strand + '1'))
				feature.qualifiers['locus_tag'] = parent
				if parent in hits:
					feature.qualifiers['note'] = hits[parent]
				features[seqname].append(feature)

with open(output_file, 'w') as out_fh:
	with open(fna_file) as in_fh:
		rec_num = 1
		for record in SeqIO.parse(in_fh, 'fasta'):
			record.annotations['molecule_type'] = 'DNA'
			record.name = "locus_%d" % rec_num
			for feature in features[record.id]:
				transl = feature.translate(record, cds = False)
				feature.qualifiers['translation'] = transl.seq.rstrip('*')
				record.features.append(feature)
			SeqIO.write(record, out_fh, 'genbank')
			rec_num += 1

