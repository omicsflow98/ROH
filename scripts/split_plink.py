#!/usr/bin/env python3

import csv
import os
import argparse
import re

parser = argparse.ArgumentParser(description="Split PLINK ROH output by chromosome")
parser.add_argument('--hom', required=True, help='Path to PLINK .hom file')
parser.add_argument('--chroms', required=True, help='Path to chromosome names file')
args = parser.parse_args()

data_file = args.hom
chrom_file = args.chroms
output_dir = "."

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

with open(data_file) as df:
    header = df.readline().strip().split()  # Read header

# Load chromosome names into a set
with open(chrom_file) as f:
    chromosomes = f.read().splitlines()
    for chrom in chromosomes:
        name = re.sub(r'[^A-Za-z0-9]', '', chrom)
        with open(data_file) as df, open(output_dir + "/" + name + ".tsv", "w") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(header)
            next(df)
            for row in df:
                row_plink = row.split()
                if row_plink[3] == chrom:
                    writer.writerow(row_plink)
