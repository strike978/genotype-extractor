#!/usr/bin/env python3
import re
import sys
import shutil
import os
from datetime import datetime

def update_ind_file(filepath):
    backup_path = filepath + '.bak'
    if os.path.exists(backup_path):
        print(f"Backup already exists: {backup_path}")
    else:
        shutil.copy2(filepath, backup_path)
        print(f"Backup: {backup_path}")
    
    output_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            indent = len(line) - len(line.lstrip())
            content = line.strip()
            
            match = re.match(r'^(\S+)\s+([MFU])\s+(.+)$', content)
            if match:
                sample_id = match.group(1)
                sex = match.group(2)
                population = match.group(3)
                
                ext = sample_id.split('.')[-1]
                new_name = f"{sample_id} {sex} {population}.{ext}"
                output_lines.append(f"{' ' * indent}{new_name}\n")
            else:
                output_lines.append(line)
    
    with open(filepath, 'w') as f:
        f.writelines(output_lines)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python update_ind.py <file.ind>")
        sys.exit(1)
    
    filepath = sys.argv[1]
    if os.path.exists(filepath):
        update_ind_file(filepath)
        print(f"Updated: {filepath}")
    else:
        print(f"File not found: {filepath}")