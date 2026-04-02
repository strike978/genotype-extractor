#!/usr/bin/env python3
# Script to generate a merged genotype CSV from 1000genomes, HGDP, and SGDP folders.
# The SNP columns are determined by the union of all SNPs present in the data files.
# Output: 113.csv

import os
import csv
from collections import defaultdict

# DATA_FOLDERS = ["1000genomes", "hgdp", "sgdp", "ggvp"]
DATA_FOLDERS = ["1000genomes", "hgdp", "sgdp"]
WORKDIR = os.path.dirname(os.path.abspath(__file__))


def get_all_sample_files():
    files = []
    for folder in DATA_FOLDERS:
        folder_path = os.path.join(WORKDIR, folder)
        if not os.path.isdir(folder_path):
            continue
        for pop in os.listdir(folder_path):
            pop_path = os.path.join(folder_path, pop)
            if not os.path.isdir(pop_path):
                continue
            for fname in os.listdir(pop_path):
                if fname.endswith(".txt"):
                    files.append((folder, pop, os.path.join(pop_path, fname)))
    return files


def collect_all_snps(files):
    snps = set()
    for _, _, fpath in files:
        with open(fpath) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    snps.add(parts[0])
    return sorted(snps)


def get_1000g_group_full_mapping():
    # Hardcoded mapping from 1000g group code to group_full, extracted from merged_genotypes.csv
    return {
        'ACB': 'African Caribbean in Barbados',
        'ASW': 'African Ancestry in Southwest US',
        'BEB': 'Bengali in Bangladesh',
        'CDX': 'Chinese Dai in Xishuangbanna, China',
        'CEU': 'Utah Residents (CEPH) with Northern and Western European Ancestry',
        'CHB': 'Han Chinese in Beijing, China',
        'CHS': 'Southern Han Chinese',
        'CLM': 'Colombian in Medellin, Colombia',
        'ESN': 'Esan in Nigeria',
        'FIN': 'Finnish in Finland',
        'GBR': 'British in England and Scotland',
        'GIH': 'Gujarati Indian in Houston, TX',
        'GWD': 'Gambian in Western Divisions in the Gambia',
        'IBS': 'Iberian Population in Spain',
        'ITU': 'Indian Telugu in the UK',
        'JPT': 'Japanese in Tokyo, Japan',
        'KHV': 'Kinh in Ho Chi Minh City, Vietnam',
        'LWK': 'Luhya in Webuye, Kenya',
        'MSL': 'Mende in Sierra Leone',
        'MXL': 'Mexican Ancestry in Los Angeles, CA',
        'PEL': 'Peruvian in Lima, Peru',
        'PJL': 'Punjabi in Lahore, Pakistan',
        'PUR': 'Puerto Rican in Puerto Rico',
        'STU': 'Sri Lankan Tamil in the UK',
        'TSI': 'Toscani in Italia',
        'YRI': 'Yoruba in Ibadan, Nigeria',
        # Add more mappings as needed from merged_genotypes.csv
    }


def get_ggvp_group_full_mapping():
    # Optionally, add mappings for GGVP populations if you want a more descriptive group_full
    # For now, just return an empty dict (can be filled in as needed)
    return {}


def parse_sample_file(fpath):
    # Returns: (individual_id, {snp: genotype})
    snp_dict = {}
    individual = None
    with open(fpath) as f:
        for line in f:
            if line.startswith("#"):
                if line.lower().startswith("#id") or line.lower().startswith("#individual") or line.lower().startswith("# sample"):
                    # Try to extract sample/individual ID from header
                    parts = line.strip().split(":")
                    if len(parts) > 1:
                        individual = parts[1].strip()
                continue
            parts = line.strip().split()
            if len(parts) >= 4:
                snp, gt = parts[0], parts[3]
                snp_dict[snp] = gt
    if not individual:
        # fallback: use filename
        individual = os.path.splitext(os.path.basename(fpath))[0]
    return individual, snp_dict


def main():
    files = get_all_sample_files()
    print(f"Found {len(files)} sample files.")
    all_snps = collect_all_snps(files)
    print(f"Total unique SNPs: {len(all_snps)}")
    group_full_map = get_1000g_group_full_mapping()

    out_path = os.path.join(WORKDIR, "genotypes.csv")
    with open(out_path, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        header = ["source", "group", "group_full", "individual"] + all_snps
        writer.writerow(header)

        ggvp_group_full_map = get_ggvp_group_full_mapping()
        for folder, pop, fpath in files:
            if folder == "1000genomes":
                source = "1000g"
                group = pop
                group_full = group_full_map.get(pop, "")
            elif folder == "hgdp":
                source = "HGDP"
                group = pop
                group_full = pop
            elif folder == "sgdp":
                source = "SGDP"
                group = pop
                group_full = pop
            # elif folder == "ggvp":
            #     source = "GGVP"
            #     group = pop
            #     group_full = ggvp_group_full_map.get(pop, pop)
            else:
                source = folder
                group = pop
                group_full = pop
            individual, snp_dict = parse_sample_file(fpath)
            row = [source, group, group_full, individual]
            for snp in all_snps:
                row.append(snp_dict.get(snp, "NA"))
            writer.writerow(row)
    print(f"Merged CSV written to {out_path}")


if __name__ == "__main__":
    main()
