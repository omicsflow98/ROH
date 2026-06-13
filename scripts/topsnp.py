#!/usr/bin/env python3

import argparse
import pandas as pd
import yaml
import sys
import os
from scipy.stats import norm

parser = argparse.ArgumentParser(description="Read counts file")
parser.add_argument('--counts', required=True, help='Path to counts file')
parser.add_argument('--type', choices=['quantile', 'pvalue'], required=True, help='type of island definition (quantile or pvalue)')
parser.add_argument('--threshold', required=True, help='threshold for defining islands (quantile or pvalue)')
parser.add_argument('--version', required=True, help='process version')
args = parser.parse_args()

countsfile = args.counts
analysis_type = args.type
threshold = float(args.threshold)
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
counts_df['z_scores'] = (counts_df[4] - mean_counts) / std_counts

counts_df['p_value'] = norm.sf(counts_df['z_scores'])
if analysis_type == 'pvalue':
    top_snps = counts_df[counts_df['p_value'] < threshold]
elif analysis_type == 'quantile':
    cutoff = counts_df['p_value'].quantile(1-threshold)
    top_snps = counts_df[counts_df['p_value'] < cutoff]

top_snps = top_snps.drop(columns=['z_scores', 'p_value'])

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