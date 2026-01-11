#!/usr/bin/env python3

import pandas as pd
import csv
import os
import argparse
import re
import sys
import yaml

parser = argparse.ArgumentParser(description="Find overlaps between ROH and genes")
parser.add_argument('--tsv', required=True, help='tsv file name')
parser.add_argument('--bed', required=True, help='bed file with gene locations')
parser.add_argument('--version', required=True, help='Nextflow process version')
args = parser.parse_args()

entry = args.tsv
bed_file = args.bed
version = args.version
output_dir = "."

with open(entry) as df:
    header = df.readline().strip().split()  # Read header
    header.append("genes")

with open(bed_file) as gb, open(entry) as rp, open(output_dir + "/" + "OL_" + entry, "w") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(header)
    next(rp)

    for row in rp:
        gb.seek(0)
        rd = csv.reader(gb, delimiter="\t")
        bed_genes = ""
        row_plink = row.split()
        row_plink[6] = int(row_plink[6])
        row_plink[7] = int(row_plink[7])
        interval_plink = pd.Interval(left=row_plink[6], right=row_plink[7])
        for row_bed in rd:
            no_special = re.sub(r'[^A-Za-z0-9]', '', row_bed[0])
            if no_special != entry.replace(".tsv", ""):
                continue
            row_bed[1] = int(row_bed[1])
            row_bed[2] = int(row_bed[2])
            interval_bed = pd.Interval(left=row_bed[1], right=row_bed[2])
            if interval_bed.overlaps(interval_plink):
                bed_genes += row_bed[3] + " "
        if bed_genes != "":
            row_plink.append(bed_genes)
        writer.writerow(row_plink)


data = {
    'Process version': {
        'overlap' : version
    },
    
    'Tool version': {
        'python' : sys.version.split()[0]
    }

}

file_path = 'versions.yml'
with open(file_path, 'w') as file:
    yaml.dump(data, file)