# extract DNA from Allen Ancient DNA Resource files

import glob
import os
import re


def process_population(selected_file):
    """Processes a selected population"""
    with open(selected_file, 'r') as f:
        data = f.readlines()

    search_term = input("Enter the population name: ")
    found_lines = [line for line in data
                   if search_term.lower() in line.split('\t')[0].lower()
                   # Filtering added here
                   and not line.split('\t')[0].startswith("Ignore_")]

    if found_lines:
        for i, line in enumerate(found_lines):
            print(f"{i+1}. {line.rstrip()}")
    else:
        print("No matches found.")
        return  # Exit the function if no matches

    # Input validation with try-except block
    while True:
        try:
            choice = int(
                input("Enter the number of the line containing the population: "))
            if choice < 1 or choice > len(found_lines):
                print("Invalid choice. Please try again.")
            else:
                break  # Valid choice received
        except ValueError:
            print("Invalid input. Please enter a number.")

    selected_line = found_lines[choice - 1]
    selected_identifier = selected_line.split(
    )[0] + '_' + selected_line.split()[1]

    # Example Operations - Modify as needed
    individuals_filename = f"{selected_identifier}_individuals.txt"
    with open(individuals_filename, 'w') as f:
        f.write(selected_line)

    plink_base_filename = os.path.splitext(selected_file)[0]
    os.system(
        f"plink --bfile {plink_base_filename} --keep {individuals_filename} --recode 23")

    try:
        os.rename("plink.txt", f"{selected_identifier}.txt")
    except FileExistsError:
        os.remove(f"{selected_identifier}.txt")
        os.rename("plink.txt", f"{selected_identifier}.txt")

    # Add other PLINK-related extensions if needed
    extensions_to_remove = ['.log', '.nosex']
    for filename in os.listdir():
        base, ext = os.path.splitext(filename)
        if ext in extensions_to_remove or base == f"{selected_identifier}_individuals":
            os.remove(filename)


# Find all .fam files in the current directory
fam_files = glob.glob("*.fam")


def flexible_sort(value):
    """Performs numeric sorting where possible, otherwise natural sorting."""
    parts = re.split(r'(\d+)', value)
    return tuple((int(part) if part.isdigit() else part) for part in parts)


fam_files.sort(key=flexible_sort)

# Display the sorted filenames
for i, filename in enumerate(fam_files):
    print(f"{i+1}. {filename}")

# Prompt the user to select a file
choice = int(input("Enter the number of the file you want to open: "))
while choice < 1 or choice > len(fam_files):
    print("Invalid choice. Please try again.")
    choice = int(input("Enter the number of the file you want to open: "))

selected_file = fam_files[choice - 1]

# Main processing loop
while True:
    process_population(selected_file)

    continue_processing = input(
        "Continue searching for populations in the same file (Y/N)? ").lower()
    if continue_processing != 'y':
        break
