
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pandas import read_csv

ref_annots = []

input_files = snakemake.input["gbk"]
data_file   = str(snakemake.input["data"])
output_file = str(snakemake.output)

data = read_csv(data_file, sep = "\t")

with open(output_file, 'w') as fd_out:
    for input_file in input_files:
        with open(input_file) as fd_in:
            for record in SeqIO.parse(fd_in, 'genbank'):
                for feature in record.features:
                    quals = feature.qualifiers
                    if feature.type == 'CDS' and 'locus_tag' in feature.qualifiers:
                        ID = feature.qualifiers['locus_tag'][0]
                        if ID:
                            row = data[data.ID == ID]
                            if 'note' not in feature.qualifiers:
                                feature.qualifiers['note'] = []
                            if not row.empty:
                                if row.Gene_Pfam.notnull().item():
                                    feature.qualifiers['note'].append('Gene by hhsearch Pfam match: %s' % row.Gene_Pfam.item())
                                    feature.qualifiers['product'] = row.Gene_Pfam.item()
                                if row.Gene_Cluster.notnull().item():
                                    feature.qualifiers['note'].append('Gene by cluster membership: %s' % row.Gene_Cluster.item())
                                    feature.qualifiers['product'] = row.Gene_Cluster.item()
                                if row['Hit.ID'].notnull().item() and row.Probab.item:
                                    feature.qualifiers['note'].append('hhsearch best hit against Pfam: %s %s (%s)' % (row['Hit.ID'].item(), row['Hit.Description'].item(), row.Probab.item()))
                SeqIO.write(record, fd_out, 'genbank')
