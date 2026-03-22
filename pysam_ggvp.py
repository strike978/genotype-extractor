#!/usr/bin/env python3
from collections import defaultdict
import csv
import os
import pysam
import urllib.request


# GGVP metadata and VCFs
GGVP_META_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/gambian_genome_variation_project/gambian_genome_variation_project.sequence.index"
GGVP_META_FILE = "gambian_genome_variation_project.sequence.index"

# Map chromosome number to GGVP VCF URL
GGVP_VCF_URLS = {
    str(i): f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/gambian_genome_variation_project/release/20200217_biallelic_SNV/ALL_GGVP.chr{i}.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz"
    for i in list(range(1, 23))
}
GGVP_VCF_URLS["X"] = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/gambian_genome_variation_project/release/20200217_biallelic_SNV/ALL_GGVP.chrX.shapeit2_integrated_snvindels_v1b_20200120.GRCh38.phased.vcf.gz"

ALLELES_CSV = "merged_snps.csv"

# Download panel file if not exists


def download_file(url, filename):
    if not os.path.isfile(filename):
        print(f"[INFO] Downloading {filename} from {url} ...")
        try:
            urllib.request.urlretrieve(url, filename)
            print(f"[INFO] Downloaded {filename}.")
        except Exception as e:
            print(f"[ERROR] Failed to download {filename}: {e}")
    else:
        print(f"[INFO] {filename} already exists. Skipping download.")

# Download VCF and index if not exists


def collect_all_snp_data(snps_by_chrom, snp_info):
    """Collect all SNP data across chromosomes and organize by sample (GGVP)"""
    print("[INFO] Collecting all SNP data from GGVP...")

    all_sample_data = defaultdict(list)
    missing_snps = []

    for idx, (chrom, positions) in enumerate(snps_by_chrom.items(), 1):
        print(
            f"\n[INFO] Processing chromosome {chrom} ({idx}/{len(snps_by_chrom)}) with {len(positions)} SNPs...")

        vcf_url = GGVP_VCF_URLS.get(chrom)
        if not vcf_url:
            print(f"[WARNING] No VCF URL for chromosome {chrom}. Skipping.")
            continue
        print(
            f"    [INFO] Fetching data for chr{chrom} from remote GGVP VCF...")
        try:
            vcf = pysam.VariantFile(vcf_url)
        except Exception as e:
            print(f"[ERROR] Could not open VCF for chr{chrom}: {e}")
            continue

        # Print VCF sample names for debugging
        vcf_samples = list(vcf.header.samples)
        print(
            f"    [DEBUG] VCF samples (first 10): {vcf_samples[:10]} ... total {len(vcf_samples)}")

        for pos_idx, pos in enumerate(positions, 1):
            print(
                f"      [PROGRESS] Fetching variant {pos_idx}/{len(positions)} at {chrom}:{pos}")
            found_variant = False
            for record in vcf.fetch(chrom, int(pos)-1, int(pos)):
                found_variant = True
                chrom_col = record.chrom
                pos_col = str(record.pos)
                id_col = record.id
                ref = record.ref
                alt = record.alts[0] if record.alts else ''

                key = (chrom_col, pos_col)
                if key in snp_info:
                    current_variant = {
                        'rsid': snp_info[key]['rsID'],
                        'chrom': snp_info[key]['chrom'],
                        'pos': snp_info[key]['pos'],
                        'ref': ref,
                        'alt': alt
                    }
                    print(
                        f"        [PROGRESS] Variant {current_variant['rsid']} at {chrom_col}:{pos_col}")

                    for sample_idx, sample in enumerate(record.samples, 1):
                        gt = record.samples[sample].get('GT')
                        if gt is None or None in gt:
                            genotype = '--'
                        else:
                            gt0, gt1 = gt
                            if gt0 == 0 and gt1 == 0:
                                genotype = ref + ref
                            elif gt0 == 1 and gt1 == 1:
                                genotype = alt + alt
                            elif (gt0 == 0 and gt1 == 1) or (gt0 == 1 and gt1 == 0):
                                genotype = ref + alt
                            else:
                                genotype = '--'

                        all_sample_data[sample].append((
                            current_variant['rsid'],
                            current_variant['chrom'],
                            current_variant['pos'],
                            genotype
                        ))
                        if sample_idx % 500 == 0:
                            print(
                                f"          [PROGRESS] {sample_idx} samples processed for variant {current_variant['rsid']}")
            if not found_variant:
                print(
                    f"      [DEBUG] No variant found at {chrom}:{pos} in VCF.")
                # Track missing SNPs by rsID and position
                for key, val in snp_info.items():
                    if key == (chrom, pos):
                        missing_snps.append(
                            {'rsid': val['rsID'], 'chrom': chrom, 'pos': pos})

        vcf.close()
        print(f"    [INFO] Chromosome {chrom} processed.")

    return all_sample_data, missing_snps


def write_sample_files(sample_data, meta_file):
    """Write individual sample files organized by population (GGVP)"""
    print("[INFO] Writing individual sample files...")

    # Load population data from GGVP metadata
    pop_info = {}
    # First, find the header line
    with open(meta_file) as f:
        for line in f:
            if line.startswith('#') and not line.startswith('##'):
                header = line.lstrip('#').strip().split('\t')
                break
        else:
            raise ValueError('No header line found in metadata file')
    # Now, parse the data lines
    sample_name_prefixes = []
    with open(meta_file) as f:
        header_found = False
        for line in f:
            if not header_found:
                if line.startswith('#') and not line.startswith('##'):
                    header_found = True
                continue
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            entry = dict(zip(header, parts))
            sample_name_full = entry.get('SAMPLE_NAME', '')
            sample_name_prefix = sample_name_full.split(
                '-')[0] if '-' in sample_name_full else sample_name_full
            pop = entry.get('POPULATION', '')
            sample_id = entry.get('SAMPLE_ID', '')
            pop_info[sample_name_prefix] = {
                'pop': pop,
                'sample_id': sample_id,
                'full_sample_name': sample_name_full
            }
            sample_name_prefixes.append(sample_name_prefix)
    print(
        f"[DEBUG] First 10 parsed sample_name_prefixes from metadata: {sample_name_prefixes[:10]}")

    # Create output directory structure for GGVP
    output_dir = "ggvp"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Print GGVP-like sample names from VCF and metadata for debugging
    vcf_samples = set(sample_data.keys())
    vcf_ggvp = sorted(
        [s for s in vcf_samples if s.startswith('SC_') or s.startswith('GM_')])
    meta_ggvp = sorted(
        [s for s in pop_info.keys() if s.startswith('SC_') or s.startswith('GM_')])
    intersection = sorted(set(vcf_ggvp) & set(meta_ggvp))
    print(f"[DEBUG] GGVP-like samples in VCF: {vcf_ggvp}")
    print(f"[DEBUG] GGVP-like samples in metadata: {meta_ggvp}")
    print(f"[DEBUG] Intersection (should be written): {intersection}")

    # Only keep samples present in both VCF and GGVP metadata
    ggvp_samples = [s for s in vcf_ggvp if s in pop_info]
    print(
        f"[INFO] Found {len(ggvp_samples)} GGVP samples present in both VCF and metadata.")

    # Group GGVP samples by population
    samples_by_pop = defaultdict(list)
    for sample_name in ggvp_samples:
        pop = pop_info[sample_name]['pop']
        samples_by_pop[pop].append(sample_name)

    # Write files
    for pop, sample_names in samples_by_pop.items():
        pop_dir = os.path.join(output_dir, pop)
        if not os.path.exists(pop_dir):
            os.makedirs(pop_dir)

        print(
            f"    [INFO] Writing {len(sample_names)} samples for population {pop}...")

        for idx, sample_name in enumerate(sample_names, 1):
            sample_file = os.path.join(pop_dir, f"{sample_name}.txt")
            with open(sample_file, 'w') as f:
                # Write sample details at the top
                pop_val = pop_info.get(sample_name, {}).get('pop', '')
                sample_id_val = pop_info.get(
                    sample_name, {}).get('sample_id', '')
                f.write(f"# Sample : {sample_name}\n")
                f.write(f"# Sample_ID : {sample_id_val}\n")
                f.write(f"# Population : {pop_val}\n")
                # Write header
                f.write("# rsid\tchromosome\tposition\tgenotype\n")

                # Sort SNPs by chromosome and position
                snp_data = sorted(sample_data[sample_name], key=lambda x: (
                    int(x[1]) if x[1].isdigit() else 999, int(x[2])))

                # Write SNP data
                for rsid, chrom, pos, genotype in snp_data:
                    f.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
            if idx % 100 == 0 or idx == len(sample_names):
                print(
                    f"      [PROGRESS] {idx}/{len(sample_names)} samples written for population {pop}")

    print(
        f"[INFO] Individual sample files written to {output_dir}/ organized by population")
    return output_dir


def main():
    print("[INFO] Script started. Creating individual sample files like 23andMe format (GGVP).")
    print("[INFO] Downloading GGVP metadata file if needed...")
    download_file(GGVP_META_URL, GGVP_META_FILE)

    print(
        f"[INFO] Reading SNP list from {ALLELES_CSV} and grouping by chromosome...")
    snps_by_chrom = defaultdict(list)
    snp_info = {}

    # Detect which columns to use for chromosome and position (chr37 or chr38)
    with open(ALLELES_CSV, 'r') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if 'chr38' in header and 'pos38' in header:
            chrom_col = 'chr38'
            pos_col = 'pos38'
            print("[INFO] Detected chr38/pos38 columns. Using GRCh38 coordinates.")
        elif 'chr37' in header and 'pos37' in header:
            chrom_col = 'chr37'
            pos_col = 'pos37'
            print("[INFO] Detected chr37/pos37 columns. Using GRCh37 coordinates.")
        else:
            raise ValueError(
                "Could not detect chr37/chr38 columns in the CSV header.")

        for row in reader:
            rsid = row['rsID']
            chrom = str(row[chrom_col])
            pos = str(row[pos_col])
            snps_by_chrom[chrom].append(pos)
            snp_info[(chrom, pos)] = {
                'rsID': rsid,
                'chrom': chrom,
                'pos': pos
            }

    print(
        f"[INFO] Loaded {len(snp_info)} SNPs across {len(snps_by_chrom)} chromosomes")

    # Collect all SNP data organized by sample
    all_sample_data, missing_snps = collect_all_snp_data(
        snps_by_chrom, snp_info)
    print(f"[INFO] Collected data for {len(all_sample_data)} samples")

    # Write individual sample files organized by population
    output_dir = write_sample_files(all_sample_data, GGVP_META_FILE)

    # Report missing SNPs at the end
    if missing_snps:
        print(
            f"\n[REPORT] {len(missing_snps)} SNPs from merged_snps.csv were not found in the GGVP VCFs:")
        for snp in missing_snps:
            print(f"  - {snp['rsid']} (chr{snp['chrom']}:{snp['pos']})")
    else:
        print("\n[REPORT] All SNPs from merged_snps.csv were found in the GGVP VCFs.")

    print(f"\n[INFO] Done! Individual sample files written to {output_dir}/")
    print(f"[INFO] Each population folder contains individual sample files in 23andMe format.")


if __name__ == "__main__":
    main()
