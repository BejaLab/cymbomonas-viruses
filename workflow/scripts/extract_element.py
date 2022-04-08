
import re
from Bio import SeqIO

input_gbk = str(snakemake.input)
output_gbk = str(snakemake.output)
flank = snakemake.params['flank']
segment = snakemake.params['segment']

coords_re = re.compile(r'(.+)_(\d+)-(\d+)[. ]*$')
def get_coords(description):
    search = coords_re.search(description)
    assert search, "Unexpected sequence name: %s" % description
    scaffold, start, end = search.groups()
    start = int(start)
    end = int(end) + flank
    if start > flank:
        start = start - flank
    else:
        start = 1
    return scaffold, start, end

found = False
with open(input_gbk) as fd:
    records = SeqIO.parse(fd, 'genbank')
    for record in records:
        scaffold, start, end = get_coords(record.description)
        if scaffold == segment['scaffold'] and start <= segment['scaffold_start'] and end >= segment['scaffold_end']:
            found = True
            break
assert found, "Segment %s not found" % segment['segment_id']

segment_start = segment['scaffold_start'] - start
segment_end   = segment['scaffold_end'] - start + 1

sub_record = record[segment_start:segment_end]
if segment['strand'] < 0:
    sub_record = sub_record.reverse_complement()

sub_record.name = segment['name']
sub_record.id = "%s_%d-%d" % (segment['scaffold'], segment['scaffold_start'], segment['scaffold_end'])
sub_record.description = sub_record.id

sub_record.annotations["molecule_type"] = "DNA"
SeqIO.write(sub_record, output_gbk, "genbank")
