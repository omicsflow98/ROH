#!/usr/bin/env python3

"""
Merge multiple versions.yml files into a single versions.yml.
Conflicts in process/tool versions raise an error.
"""

import yaml
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="Merge Nextflow versions.yml files")
parser.add_argument("--files", required=True,  nargs="+", help="Paths to versions.yml files to merge")
parser.add_argument("-o", "--output", default="version_info.yml", help="Output merged YAML file (default: versions.yml)")
args = parser.parse_args()

files = args.files
output = args.output

merged = defaultdict(dict)

for vf in files:
    with open(vf) as f:
        data = yaml.safe_load(f) or {}

    for section, entries in data.items():
        merged.setdefault(section, {})

        for key, value in entries.items():
            if key in merged[section]:
                if merged[section][key] != value:
                    # Append new value (avoid duplicates)
                    existing = merged[section][key]

                    # Turn existing into a list if it's not already
                    if not isinstance(existing, list):
                        existing = [v.strip() for v in str(existing).split(",")]

                    # Add new value if it's not already present
                    if str(value) not in existing:
                        existing.append(str(value))

                    # Store back as comma-separated string
                    merged[section][key] = ", ".join(existing)
            else:
                merged[section][key] = value

with open(output, "w") as out:
    yaml.dump(dict(merged), out, sort_keys=False, default_flow_style=False)
