#!/usr/bin/env python3

import argparse
import sys
import yaml
import pyranges as pr
import pandas as pd

parser = argparse.ArgumentParser(description="Find overlaps between ROH and genes (pure PyRanges)")
parser.add_argument('--tsv', required=True, help='tsv file name')
parser.add_argument('--bed', required=True, help='bed file with gene locations')
parser.add_argument('--version', required=True, help='Nextflow process version')
args = parser.parse_args()

entry = args.tsv
bed_file = args.bed
version = args.version
output_dir = "."

chrom_name = entry.replace(".tsv", "")

# -------------------------
# Load ROH as PyRanges
# -------------------------
roh_df = pd.read_csv(entry, sep="\t")

# Add required columns
roh_df["Chromosome"] = chrom_name
roh_df["Start"] = roh_df.iloc[:, 6].astype(int)
roh_df["End"] = roh_df.iloc[:, 7].astype(int)
roh_df["__idx__"] = roh_df.index

roh = pr.PyRanges(roh_df[["Chromosome", "Start", "End", "__idx__"]])

# -------------------------
# Load BED as PyRanges
# -------------------------
genes = pr.read_bed(bed_file)

# Clean chromosome names to match your logic
genes_df = genes.df
genes_df["Chromosome"] = genes_df["Chromosome"].astype(str).str.replace(r'[^A-Za-z0-9]', '', regex=True)

# Filter to relevant chromosome
genes_df = genes_df[genes_df["Chromosome"] == chrom_name]

genes = pr.PyRanges(genes_df)

# -------------------------
# Find overlaps
# -------------------------
overlaps = roh.join(genes)

# -------------------------
# Aggregate genes per ROH
# -------------------------
if overlaps.df.empty:
    gene_map = {}
else:
    gene_map = (
        overlaps.df
        .groupby("__idx__")["Name"]  # BED gene column = Name in PyRanges
        .apply(lambda x: " ".join(sorted(set(x))))
        .to_dict()
    )

# -------------------------
# Write output manually
# -------------------------
with open(entry) as rp, open(output_dir + "/" + "OL_" + entry, "w") as out:
    header = rp.readline().strip().split("\t")
    header.append("genes")
    out.write("\t".join(header) + "\n")

    for i, line in enumerate(rp):
        row = line.strip().split("\t")
        genes = gene_map.get(i, "")
        row.append(genes)
        out.write("\t".join(row) + "\n")

# -------------------------
# Write version info
# -------------------------
data = {
    'Process version': {
        'overlap': version
    },
    'Tool version': {
        'python': sys.version.split()[0]
    }
}

with open('versions.yml', 'w') as file:
    yaml.dump(data, file)