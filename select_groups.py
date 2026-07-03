import sys
import os

FAM_FILE = "v66.HO.aadr.PUB.fam"
KEEP_FILE = "keep.txt"


def read_fam_groups(fam_path):
    groups = {}
    with open(fam_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            group, sample = parts[0], parts[1]
            if group not in groups:
                groups[group] = []
            groups[group].append(sample)
    return groups


def search_groups(groups, query):
    query = query.lower()
    return [g for g in groups if query in g.lower()]


def main():
    fam_path = FAM_FILE
    if not os.path.exists(fam_path):
        print(f"FAM file not found: {fam_path}")
        sys.exit(1)
    groups = read_fam_groups(fam_path)
    selected = set()
    while True:
        query = input(
            "Search population group (or leave blank to finish): ").strip()
        if not query:
            break
        matches = search_groups(groups, query)
        if not matches:
            print("No groups found. Try again.")
            continue
        print("Matching groups:")
        for idx, g in enumerate(matches):
            print(f"  {idx+1}. {g}")
        sel = input(
            "Select group numbers to add (comma separated, or blank to skip): ").strip()
        if sel:
            try:
                nums = [int(s)-1 for s in sel.split(",")
                        if s.strip().isdigit()]
                for n in nums:
                    if 0 <= n < len(matches):
                        selected.add(matches[n])
                if selected:
                    print("Selected populations so far:")
                    for pop in sorted(selected):
                        print(f"  - {pop}")
            except Exception:
                print("Invalid selection. Try again.")
    if not selected:
        print("No groups selected. Exiting.")
        return
    # Collect all (group, sample) pairs for selected groups
    keep_entries = []
    with open(fam_path, "r") as fam:
        for line in fam:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            group, sample = parts[0], parts[1]
            if group in selected:
                keep_entries.append((group, sample))
    # Sort by group, then by sample format (e.g., .DG, .HO), then by sample name

    def sort_key(entry):
        group, sample = entry
        # Extract format (e.g., .DG, .HO, etc.)
        if '.' in sample:
            fmt = sample.split('.')[-1]
        else:
            fmt = ''
        return (group, fmt, sample)
    keep_entries.sort(key=sort_key)
    with open(KEEP_FILE, "w") as out:
        for group, sample in keep_entries:
            out.write(f"{group}\t{sample}\n")
    print(
        f"Wrote samples for {len(selected)} groups to {KEEP_FILE} in PLINK keep format.")


if __name__ == "__main__":
    main()
