
from collections import defaultdict
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Record import HSP
from BCBio import GFF
from urllib.parse import unquote
from io import StringIO

xml_file   = snakemake.input['xml']
fasta_file = snakemake.input['fasta']
query_file = snakemake.input['query']
gff_file   = str(snakemake.output)

query_gap   = snakemake.params['query_gap']
sbjct_gap_1 = snakemake.params['sbjct_gap_1']
sbjct_gap_2 = snakemake.params['sbjct_gap_2']
chain_len   = snakemake.params['chain_len']
min_neighbor_len = snakemake.params['min_neighbor_len']

class HSP2(HSP):
    def append_next(self, hsp_next):
        if not hasattr(self, 'next'):
            self.next = []
        self.next.append(hsp_next)
    
    def join(self, chain_id):
        if hasattr(self, 'chain_id'):
            if self.chain_id != chain_id:
                return []
        else:
            self.chain_id = chain_id

        next_hsp = None
        if hasattr(self, 'next'):
            for hsp in self.next:
                if not hasattr(hsp, 'chain_id') or hsp.chain_id == chain_id:
                    if next_hsp is None or next_hsp.align_length < hsp.align_length:
                        next_hsp = hsp
        if next_hsp is None:
            return [ self ]
        else:
            return [ self ] + next_hsp.join(chain_id)

    def get_strand(self):
        return self.sbjct_start < self.sbjct_end

    def get_start(self):
        if self.sbjct_start < self.sbjct_end:
            return self.sbjct_start
        else:
            return self.sbjct_end

    def get_end(self):
        if self.sbjct_start < self.sbjct_end:
            return self.sbjct_end
        else:
            return self.sbjct_start

    def is_neighbor(self, hsp2):
        global query_gap, sbjct_gap_1, sbjct_gap_2
        if self.align_length < min_neighbor_len or hsp2.align_length < min_neighbor_len or hsp2.query_start < self.query_end or hsp2.query_start > self.query_end + query_gap:
            return False
        strand = self.get_strand()
        if strand != hsp2.get_strand():
            return False
        s1 = self.sbjct_start
        e1 = self.sbjct_end
        s2 = hsp2.sbjct_start
        e2 = hsp2.sbjct_end
        return strand == True and s2 > e1 and s2 < e1 + sbjct_gap_1 or strand == False and s1 > e2 and s1 < e2 + sbjct_gap_1

    def is_close(self, hsp2):
        global query_gap, sbjct_gap_1, sbjct_gap_2
        if not self.is_neighbor(hsp2):
            return False
        strand = self.get_strand()
        s1 = self.sbjct_start
        e1 = self.sbjct_end
        s2 = hsp2.sbjct_start
        e2 = hsp2.sbjct_end
        return strand == True and s2 < e1 + sbjct_gap_2 or strand == False and s1 < e2 + sbjct_gap_2

    def merge(self, hsp2):
        strand = self.get_strand()
        if strand == True:
            self.sbjct_start = min(self.sbjct_start, hsp2.sbjct_start)
            self.sbjct_end   = max(self.sbjct_end,   hsp2.sbjct_end)
        else:
            self.sbjct_start = max(self.sbjct_start, hsp2.sbjct_start)
            self.sbjct_end   = min(self.sbjct_end,   hsp2.sbjct_end)
        self.query_start = min(self.query_start, hsp2.query_start)
        self.query_end   = max(self.query_end,   hsp2.query_end)
        self.identities   += hsp2.identities
        self.align_length += hsp2.align_length
        self.match        += hsp2.match
        self.positives    += hsp2.positives
        self.score        += hsp2.score
        if hasattr(hsp2, 'next'):
            self.next = hsp2.next 

class Chain():
    def __init__(self, first_hsp, hit_def, query, chain_id):
        self.hsps = first_hsp.join(chain_id)
        self.id = chain_id
        self.query = query
        self.hit_def = hit_def
        self.strand = first_hsp.get_strand()
        self.align_length = 0
        for hsp in self.hsps:
            self.align_length += hsp.align_length

    def collapse(self):
        hsp = self.hsps.pop(0)
        new_hsps = [ hsp ]
        self.start = hsp.get_start()
        self.end   = hsp.get_end()
        for hsp_next in self.hsps:
            self.start = min(self.start, hsp_next.get_start())
            self.end   = max(self.end,   hsp_next.get_end())
            if hsp.is_close(hsp_next):
                hsp.merge(hsp_next)
            else:
                new_hsps.append(hsp_next)
                hsp = hsp_next
        self.hsps = new_hsps

    def to_feature(self):
        strand = 1 if self.strand else -1
        features = []
        i = 0
        for hsp in self.hsps:
            qualifiers = { 'source': "blastn", 'score': self.align_length / len(self.query.seq), 'ID': "%s-%s-%s-%s" % (self.hit_def, self.strand, self.query.id, self.id) }
            feature = SeqFeature(FeatureLocation(hsp.get_start() - 1, hsp.get_end()), type = "fragment", strand = strand, qualifiers = qualifiers)
            features.append(feature)
            i += 1
        return features

all_hsps = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

sequences = {}
with open(fasta_file) as fd:
    records = SeqIO.parse(fd, 'fasta')
    for record in records:
        sequences[record.id] = record

query_records = {}
with open(query_file) as fd:
    records = SeqIO.parse(fd, 'fasta')
    for record in records:
        query_records[record.id] = record

with open(xml_file) as fd:
    records = NCBIXML.parse(fd)
    for record in records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                strand = hsp.strand[1]
                hsp.__class__ = HSP2
                all_hsps[alignment.hit_def][strand][record.query].append(hsp)

records = []
with StringIO() as gff:
    for hit_def, strands in all_hsps.items():
        chain_id = 0
        record = sequences[hit_def] 
        features = []
        for strand, queries in strands.items():
            for query, hsps in queries.items():
                hsps = sorted(hsps, key = lambda hsp: hsp.query_start) 
                for hsp1 in hsps:
                    for hsp2 in hsps:
                        if hsp1.is_neighbor(hsp2):
                            hsp1.append_next(hsp2)
                chains = []
                for hsp in hsps:
                    chain = Chain(hsp, hit_def, query_records[query], chain_id)
                    chain_id += 1
                    if chain.align_length >= chain_len:
                        chain.collapse()
                        chains.append(chain)
                for chain in chains:
                    record.features += chain.to_feature()
        record.features = sorted(record.features, key = lambda feature: feature.location.start)
        GFF.write([record], gff)
    with open(gff_file, 'w') as fd:
        for line in gff.getvalue().split('\n'):
            if line and not line.startswith('#'):
                fd.write(line + '\n')
