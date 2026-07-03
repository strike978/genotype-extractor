import subprocess
import os
import shutil
import sys
from pathlib import Path

PLINK = "plink"

# Input files
BFILE = "v66.HO.aadr.PUB"
KEEP = "keep.txt"

# Output prefixes
SUBSET = "subset"
FILTERED = "subset.filtered"
MAF = "subset.filtered.maf"
PRUNED = "pruned"
FINAL = "subset.final"

# Keep only final files? If True, intermediate files are deleted after the pipeline completes.
CLEANUP_INTERMEDIATE = True

# Filtering parameters
GENO = "0.05"  # SNP missingness
MIND = "0.05"  # Individual missingness
MAF_THR = "0.05"  # Minor allele frequency
LD_WINDOW = "50"
LD_STEP = "5"
LD_R2 = "0.2"


def run(cmd_args):
    # cmd_args: list of command arguments (safer than shell=True)
    print(f"\n[Running] {' '.join(cmd_args)}")
    subprocess.run(cmd_args, check=True)


def main():
    # Verify PLINK is available
    if shutil.which(PLINK) is None:
        print(
            f"Error: '{PLINK}' not found on PATH. Please install PLINK or adjust PLINK variable.")
        sys.exit(1)
    # 1. Subset samples
    run([PLINK, "--bfile", BFILE, "--keep", KEEP, "--make-bed", "--out", SUBSET])
    # 2. Filter missingness
    run([PLINK, "--bfile", SUBSET, "--geno", GENO,
        "--mind", MIND, "--make-bed", "--out", FILTERED])
    # 3. Filter by MAF
    run([PLINK, "--bfile", FILTERED, "--maf", MAF_THR, "--make-bed", "--out", MAF])
    # 4. LD prune
    run([PLINK, "--bfile", MAF, "--indep-pairwise",
        LD_WINDOW, LD_STEP, LD_R2, "--out", PRUNED])
    run([PLINK, "--bfile", MAF, "--extract",
        f"{PRUNED}.prune.in", "--make-bed", "--out", FINAL])
    if CLEANUP_INTERMEDIATE:
        intermediate_prefixes = [SUBSET, FILTERED, MAF]
        for prefix in intermediate_prefixes:
            for ext in [".bed", ".bim", ".fam", ".log", ".hh", ".nosex", ".irem"]:
                p = Path(f"{prefix}{ext}")
                if p.exists():
                    p.unlink()
        # Remove non-essential final files but keep final .bed/.bim/.fam
        for ext in [".log", ".hh", ".nosex"]:
            p = Path(f"{FINAL}{ext}")
            if p.exists():
                p.unlink()
        # Prune helper files
        for pth in [f"{PRUNED}.prune.in", f"{PRUNED}.prune.log", f"{PRUNED}.prune.out", f"{PRUNED}.hh", f"{PRUNED}.log.txt", f"{PRUNED}.log"]:
            p = Path(pth)
            if p.exists():
                p.unlink()
    print("\nPLINK filtering pipeline complete. Check output files for results.")


if __name__ == "__main__":
    main()
