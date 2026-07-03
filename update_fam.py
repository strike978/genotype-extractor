#!/usr/bin/env python3
import re
import sys
import shutil
import os

def update_fam(fam_path, ind_path):
    backup_path = fam_path + '.bak'
    if os.path.exists(backup_path):
        print(f"Backup already exists: {backup_path}")
    else:
        shutil.copy2(fam_path, backup_path)
        print(f"Backup: {backup_path}")
    
    sample_map = {}
    with open(ind_path, 'r') as f:
        for line in f:
            content = line.strip()
            match = re.match(r'^(\S+)\s+([MFU])\s+(.+)$', content)
            if match:
                sample_id = match.group(1)
                population = match.group(3)
                sample_map[sample_id] = population
    
    output_lines = []
    with open(fam_path, 'r') as f:
        for line in f:
            indent = len(line) - len(line.lstrip())
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] in sample_map:
                parts[0] = sample_map[parts[1]]
            output_lines.append(' ' * indent + ' '.join(parts) + '\n')
    
    with open(fam_path, 'w') as f:
        f.writelines(output_lines)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python update_fam.py <fam_file> <ind_file>")
        sys.exit(1)
    
    fam_path = sys.argv[1]
    ind_path = sys.argv[2]
    
    if os.path.exists(fam_path) and os.path.exists(ind_path):
        update_fam(fam_path, ind_path)
        print(f"Updated: {fam_path}")
    else:
        print("File not found")