from igraph import Graph
from pybedtools import BedTool
from Bio import SeqIO
from collections import defaultdict
from os.path import getsize
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import re

gff_input_file  = str(snakemake.input['gff'])
fna_input_file  = str(snakemake.input['fna'])
gff_output_file = str(snakemake.output['gff'])
fna_output_file  = str(snakemake.output['fna'])

def get_chain_groups(chain_edges):
    g = Graph.TupleList(chain_edges)
    membership = g.components().membership
    chain_groups = defaultdict(list)
    for i, v in enumerate(g.vs):
        chain = v.attributes()['name']
        group = membership[i]
        chain_groups[group].append(chain)
    return chain_groups

id_re = re.compile('ID=([^ ;]+)')
with open(gff_output_file, 'w') as gff_out:
    with open(fna_output_file, 'w') as fna_out:
        if getsize(gff_input_file) > 10:
            records = {}
            with open(fna_input_file) as fna_in:
                fna = SeqIO.parse(fna_in, 'fasta')
                for record in fna:
                    records[record.id] = record
            chain_scores = {}
            chain_edges = set()
            intervals = {}
            chain_intervals = defaultdict(list)
            gff = BedTool(gff_input_file)
            bed = gff.merge(s = True, c = '7,6,9', o = 'distinct,collapse,collapse')
            for line in bed:
                sbjct, offset, end, strand, scores, chains = line.fields
                offset = int(offset)
                end = int(end)
                scores = scores.split(',')
                chains = chains.split(',')
                chain_len = len(chains)
                interval = sbjct, offset, end, strand
                for i, chain in enumerate(chains):
                    chain_scores[chain] = float(scores[i])
                    chain_intervals[chain].append(interval)
                    for j in range(i + 1, chain_len):
                        chain_edges.add((chains[i], chains[j]))
            chain_groups = get_chain_groups(chain_edges)
            for group, chains in chain_groups.items():
                highest_score = 0
                highest_chain = ''
                intervals = set()
                for chain in chains:
                    if chain_scores[chain] > highest_score:
                        highest_score = chain_scores[chain]
                        highest_chain = chain
                    for interval in chain_intervals[chain]:
                        intervals.add(interval)
                intervals = sorted(intervals, key = lambda interval: interval[1])
                locs = []
                desc = []
                for interval in intervals:
                    sbjct, offset, end, strand = interval
                    gff_line = [ sbjct, 'blastn', 'fragment', str(offset + 1), str(end), str(highest_score), strand, '.', highest_chain ]
                    gff_out.write('\t'.join(gff_line) + '\n')
                    desc.append('%d-%d' % (offset, end))
                    locs.append(FeatureLocation(offset + 1, end, -1 if strand == '-' else 1))
                search = id_re.search(highest_chain)
                strand = intervals[0][3]
                if len(locs) > 1:
                    if strand == '-': locs.reverse()
                    feature = SeqFeature(CompoundLocation(locs))
                else:
                    feature = SeqFeature(locs[0])
                seq = feature.extract(records[sbjct])
                seq.id = seq.name = search.group(1)
                seq.description = ''
                # seq.description = 'scaffold:%s strand:%s score:%s intervals:[%s]' % (sbjct, strand, highest_score, ','.join(desc))
                SeqIO.write(seq, fna_out, 'fasta')
