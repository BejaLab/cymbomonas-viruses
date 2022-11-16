
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from collections import defaultdict

ref_annots = []

segments_gff_files = str(snakemake.input["gff"])
proteins_gff_files = str(snakemake.input["fna"])


blastn_file = str(snakemake.input["blastn"])
pfam_file = str(snakemake.input["pfam"])
hmm_dbs = snakemake.input["hmm_dbs"]
hmmscan_files = snakemake.input["hmmscan"]
markers_files = snakemake.input["markers"]

output_file = str(snakemake.output)

evalue_threshold = snakemake.params["evalue"]
coding_feature = snakemake.params["coding_feature"]

min_repeat_id = snakemake.params["min_repeat_id"]
min_repeat_gap = snakemake.params["min_repeat_gap"]
min_repeat_len = snakemake.params["min_repeat_len"]

aliases = {
    "MCP": "Major capsid protein (MCP)",
    "mCP": "Minor capsid protein (mCP)",
    "uncharacterized protein": "-",
    "DNA polymerase - helicase": "DNA polymerase-helicase"
}

features = defaultdict(list)

with open(blastn_file) as fh:
    for line in fh:
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.rstrip().split()
        qstart = int(qstart)
        qend   = int(qend)
        sstart = int(sstart)
        send   = int(send)
        if float(pident) > min_repeat_id and qseqid == sseqid and int(length) >= min_repeat_len and qstart < qend and sstart > send and sstart - qend >= min_repeat_gap:
            loc1 = FeatureLocation(qstart - 1, qend, strand = 1)
            loc2 = FeatureLocation(send - 1, sstart, strand = -1)
            feature = SeqFeature(CompoundLocation([loc1, loc2]), type = 'repeat_region')
            feature.qualifiers['rpt_type'] = 'inverted'
            feature.qualifiers['inference'] = 'ab initio prediction:blastn'
            feature.qualifiers['note'] = 'pident: %s%%' % pident
            features[qseqid].append(feature)

descs = {}
for fname in hmm_dbs:
    with open(fname) as fh:
        NAME = ACC = None
        for line in fh:
            if line.startswith('ACC'):
                key, ACC = line.rstrip().split()
            if line.startswith('NAME'):
                key, NAME = line.rstrip().split()
            if line.startswith('DESC'):
                key, DESC = line.rstrip().split(maxsplit = 1)
                if ACC: descs[ACC] = DESC
                if NAME: descs[NAME] = DESC
                NAME = ACC = None

hits = defaultdict(list)

with open(pfam_file) as fh:
    for line in fh:
        if not line.startswith('#') and line != '\n':
            seq_id, aln_start, aln_end, env_start, env_end, hmm_acc, hmm_name, type, hmm_start, hmm_end, hmm_len, bit_score, E_value, significance, clan = line.rstrip().split()
            hit = "%s: %s (%s)" % (hmm_acc, descs[hmm_acc], E_value)
            hits[seq_id].append(hit)

for fname in markers_files:
    with open(fname) as fh:
        for line in fh:
            if not line.startswith('#'):
                t_name, t_acc, q_name, q_acc, E_value, score, bias, best_E_value, best_score, best_bias, exp, reg, clu, ov, env, dom, rep, inc, desc = line.rstrip().split(maxsplit = 18)
                if float(E_value) <= evalue_threshold:
                    hit = "Marker %s: %s (%s)" % (q_name, descs[q_name], E_value)
                    hits[t_name].append(hit)

for fname in hmmscan_files:
    with open(fname) as fh:
        for line in fh:
            if not line.startswith('#'):
                t_name, t_acc, q_name, q_acc, E_value, score, bias, best_E_value, best_score, best_bias, exp, reg, clu, ov, env, dom, rep, inc, desc = line.rstrip().split(maxsplit = 18)
                if desc in aliases:
                    desc = aliases[desc]
                definition = ''
                if t_name != desc and desc != '-':
                    definition = ' ' + desc if t_name != desc and desc != '-' else ''
                if float(E_value) <= evalue_threshold:
                    hit = "%s:%s (%s)" % (t_name, definition, E_value)
                    hits[q_name].append(hit)

with open(gff_file) as fh:
    for line in fh:
        if not line.startswith('#') and line != '\n':
            seqname, source, feature_type, start, end, score, strand, frame, attributes = line.rstrip().split('\t')
            if feature_type == coding_feature:
                atts = {}
                for att in attributes.split(';'):
                    key, value = att.split('=')
                    atts[key] = value
                if 'Parent' in atts:
                    parent = atts['Parent']
                else:
                    parent = atts['ID']
                loc = FeatureLocation(int(start) - 1, int(end))
                feature = SeqFeature(loc, type = 'CDS', strand = int(strand + '1'))
                feature.qualifiers['locus_tag'] = parent
                feature.qualifiers['inference'] = "ab initio prediction:" + source.replace("_v", ":")
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
                if feature.type == 'CDS':
                    transl = feature.translate(record, cds = False)
                    feature.qualifiers['translation'] = transl.seq.rstrip('*')
                record.features.append(feature)
            SeqIO.write(record, out_fh, 'genbank')
            rec_num += 1

