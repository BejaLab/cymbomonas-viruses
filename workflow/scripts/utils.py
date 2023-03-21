from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import itertools
import re

class GFF:
    def __init__(self, sbjct, score, strand, id):
        self.id = id
        self.sbjct = sbjct
        self.score = score
        self.strand = strand
        self.locs = []
        self.len = 0

    def add_interval(self, start, end):
        self.locs.append((start, end))
        self.len += end - start + 1

    @classmethod
    def range_diff(cls, r1, r2):
        s1, e1 = r1
        s2, e2 = r2
        endpoints = sorted((s1, s2 - 1, e1, e2 + 1))
        result = []
        if endpoints[0] == s1 and endpoints[1] != s1:
            result.append((endpoints[0], endpoints[1]))
        if endpoints[3] == e1 and endpoints[2] != e1:
            result.append((endpoints[2], endpoints[3]))
        return result

    @classmethod
    def loc_diff(cls, locs1, locs2):
        sub_locs = locs1
        for loc2 in locs2:
            sub_locs = list(itertools.chain(*[cls.range_diff(loc1, loc2) for loc1 in sub_locs]))
        return sub_locs

    @classmethod
    def locs_len(cls, locs):
        return sum([ end - start + 1 for start, end in locs ])

    def contains(self, other):
        sub_locs = self.loc_diff(other.locs, self.locs)
        sub_len = self.locs_len(sub_locs)
        return sub_len < 10

    def overlaps(self, other):
        sub_locs = self.loc_diff(self.locs, other.locs)
        return self.len - self.locs_len(sub_locs) > 10

    def add(self, other):
        intervals = sorted(other.locs + self.locs)
        stack = [ intervals.pop(0) ]
        for start, end in intervals:
            start0, end0 = stack[-1]
            if start0 <= start <= end0:
                stack[-1] = start0, max(end, end0)
            else:
                stack.append((start, end))
        self.locs = stack
        self.len = sum([ end - start + 1 for start, end in stack ])

    def subtract(self, other):
        self.locs = self.loc_diff(self.locs, other.locs)
        new_locs = []
        old_len = self.len
        self.len = 0
        for start, end in self.locs:
            loc_len = end - start + 1
            if loc_len >= 10:
                new_locs.append((start, end))
                self.len += loc_len
        self.locs = new_locs
        self.score = self.score / old_len * self.len

    def get_feature(self):
        feat_locs = []
        for start, end in self.locs:
            feat_locs.append(FeatureLocation(start, end, -1 if self.strand == '-' else 1))
        if len(feat_locs) > 1:
            if self.strand == '-': feat_locs.reverse()
            loc = CompoundLocation(feat_locs)
        else:
            loc = feat_locs[0]
        return SeqFeature(loc)

    def get_gff(self):
        lines = []
        for start, end in self.locs:
            line = [ self.sbjct, 'blastn', 'fragment', str(start), str(end), str(self.score), self.strand, '.', 'ID=%s;' % self.id ]
            lines.append('\t'.join(line))
        return '\n'.join(lines) + '\n'

    @classmethod
    def parse_gff(cls, fh):
        records = {}
        prev_sbjct = ''
        id_re = ''
        for line in fh:
            if line and line != '\n' and not line.startswith('#'):
                sbjct, source, type, start, end, score, strand, frame, comment = line.strip().split('\t')
                search = re.search('ID=([^ ;]+)', comment)
                assert search, 'Unexpected comment field format: %s' % comment
                id = search.group(1)
                if prev_sbjct and prev_sbjct != sbjct:
                    yield prev_sbjct, records.values()
                    records = {}
                if id not in records:
                    score = float(score) if score != '.' else -1
                    records[id] = GFF(sbjct, float(score), strand, id)
                records[id].add_interval(int(start), int(end))
                prev_sbjct = sbjct
        yield prev_sbjct, records.values()
