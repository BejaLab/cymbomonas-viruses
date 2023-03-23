from Bio import SeqIO
from utils import GFF

gff_input_file  = str(snakemake.input['gff'])
fna_input_file  = str(snakemake.input['fna'])
output_file = str(snakemake.output)

parse = SeqIO.parse(fna_input_file, 'fasta')
fasta = SeqIO.to_dict(parse)

with open(gff_input_file) as gff:
    with open(output_file, 'w') as fna:
        for sbjct, records in GFF.parse_gff(gff):
            for record in records:
                feature = record.get_feature()
                seq = feature.extract(fasta[sbjct])
                seq.id = seq.name = record.id
                seq.description = 'score:%f' % record.score
                SeqIO.write(seq, fna, 'fasta')
