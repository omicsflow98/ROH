#!/usr/bin/env python3

import argparse
import pandas as pd
import yaml
import sys
import os

parser = argparse.ArgumentParser(description="Read counts file")
parser.add_argument('--counts', required=True, help='Path to counts file')
parser.add_argument('--version', required=True, help='process version')
args = parser.parse_args()

countsfile = args.counts
version = args.version

if os.path.getsize(countsfile) == 0:
    print("No SNPs found, skipping island detection")
    data = {
        'Process version': {
            'topsnp' : version
        },
            
        'Tool version': {
            'python' : sys.version.split()[0]
        }
    }
    file_path = 'versions.yml'
    with open(file_path, 'w') as file:
        yaml.dump(data, file)
    open("topsnps.tsv", "w").close()
    sys.exit(0)
    

counts_df = pd.read_csv(countsfile, sep='\t', header=None)
mean_counts = counts_df[4].mean()
std_counts = counts_df[4].std()
z_scores = (counts_df[4] - mean_counts) / std_counts

threshold = z_scores.quantile(0.999)
top_snps = counts_df[z_scores > threshold]

top_snps.to_csv('topsnps.tsv', sep="\t", index=False, header=False)

data = {
    'Process version': {
        'topsnp' : version
    },
    
    'Tool version': {
        'python' : sys.version.split()[0]
    }

}

file_path = 'versions.yml'
with open(file_path, 'w') as file:
    yaml.dump(data, file)