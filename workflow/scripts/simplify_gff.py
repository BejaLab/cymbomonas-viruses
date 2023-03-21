from copy import deepcopy
from utils import GFF

input_file  = str(snakemake.input)
output_file = str(snakemake.output)

with open(input_file) as gff_fh:
    with open(output_file, 'w') as gff_out:
        for sbjct, records in GFF.parse_gff(gff_fh):
            records = sorted(records, key = lambda r: -r.len)
            contained = {}
            num_recs = len(records)
            # note records contained in other records
            for i in range(num_recs - 1):
                for j in range(i + 1, num_recs):
                    if j != i and j not in contained and records[i].contains(records[j]):
                        contained[j] = i

            # take only those records that are not contained in others
            uniq_records = [ records[i] for i in range(num_recs) if i not in contained ]
            uniq_records = sorted(uniq_records, key = lambda r: -r.len)
            incl_records = [ uniq_records.pop(0) ]
            mega_record = deepcopy(incl_records[0])
            to_toss = len(uniq_records)

            # add non-overlapping records
            while uniq_records and to_toss > -1:
                record = uniq_records.pop(0)
                if not mega_record.contains(record):
                    if mega_record.overlaps(record):
                        if not uniq_records:
                            to_toss = -1
                        else:
                            to_toss += -1
                        # print('Returned %s %d' % (record.id, to_toss))
                        uniq_records.append(record)
                    else:
                        mega_record.add(record)
                        incl_records.append(record)
                        # print('Included %s' % record.id)

            # add the rest of the records that are not entirely covered already
            uniq_records = sorted(uniq_records, key = lambda r: -r.len)

            for record in uniq_records:
                if not mega_record.contains(record):
                    record.subtract(mega_record)
                    if record.len > 10:
                        mega_record.add(record)
                        incl_records.append(record)

            incl_records = sorted(incl_records, key = lambda r: r.locs[0][0])

            for record in incl_records:
                gff_lines = record.get_gff()
                gff_out.write(gff_lines)
