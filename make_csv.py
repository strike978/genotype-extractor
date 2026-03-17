#!/usr/bin/env python3
# Script to generate a merged genotype CSV from 1000genomes, HGDP, and SGDP folders.
# The SNP columns are determined by the union of all SNPs present in the data files.
# Output: 113.csv

import os
import csv
from collections import defaultdict

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
    # Build mapping from group code to group_full using merged_genotypes.csv
    mapping = {}
    merged_path = os.path.join(WORKDIR, "merged_genotypes.csv")
    if not os.path.exists(merged_path):
        return mapping
    with open(merged_path) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if row[0] == "1000g":
                code = row[1]
                full = row[2]
                if code not in mapping:
                    mapping[code] = full
    return mapping


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

    out_path = os.path.join(WORKDIR, "113.csv")
    with open(out_path, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        header = ["source", "group", "group_full", "individual"] + all_snps
        writer.writerow(header)

        for folder, pop, fpath in files:
            source = ("1000g" if folder == "1000genomes" else (
                "HGDP" if folder == "hgdp" else "SGDP"))
            group = pop
            if folder == "1000genomes":
                group_full = group_full_map.get(pop, f"{pop}")
            else:
                group_full = pop  # Could be improved with metadata if needed
            individual, snp_dict = parse_sample_file(fpath)
            row = [source, group, group_full, individual]
            for snp in all_snps:
                row.append(snp_dict.get(snp, "NA"))
            writer.writerow(row)
    print(f"Merged CSV written to {out_path}")


if __name__ == "__main__":
    main()
