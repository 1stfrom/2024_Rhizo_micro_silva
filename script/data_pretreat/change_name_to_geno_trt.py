import pandas as pd
import os
import glob
from collections import defaultdict

# Path to your CSV file and the directory containing the gene files
csv_file_path = '/Users/nathanma/Documents/PHDLife/mike_proj/data/detail_info.csv'
gene_files_directory = '/Users/nathanma/Documents/PHDLife/mike_proj/data/largedata/R'
df = pd.read_csv('/Users/nathanma/Documents/PHDLife/mike_proj/data/detail_info.csv')

# Clean potential whitespace issues
df[df.columns[0]] = df[df.columns[0]].str.strip()
df[df.columns[1]] = df[df.columns[1]].str.strip()
df[df.columns[2]] = df[df.columns[2]].str.strip()

# Initialize a counter dictionary to track occurrences of A_B combinations
counter_dict = defaultdict(int)

# Function to generate unique file name based on A_B with handling duplicates
def generate_unique_name(row):
    base_name = f"{row[1]}_{row[2]}"
    # Increment the counter for this base_name and append it if > 0
    counter_dict[base_name] += 1
    counter = counter_dict[base_name]
    if counter > 1:
        unique_name = f"{row[1]}_{counter-1}_{row[2]}"
    else:
        unique_name = base_name
    return row[0], unique_name

# Generate mapping
mapping = dict(df.apply(lambda row: generate_unique_name(row), axis=1))

# Function to rename files based on the mapping and file type
def rename_files(file_type):
    for file_path in glob.glob(os.path.join(gene_files_directory, f'SRR*_{file_type}.fastq.gz')):
        file_name = os.path.basename(file_path)
        prefix = file_name.split('_')[0]
        if prefix in mapping:
            new_file_name = f"{mapping[prefix]}_{file_type}.fastq.gz"
            new_file_path = os.path.join(gene_files_directory, new_file_name)
            os.rename(file_path, new_file_path)
            print(f"Renamed {file_path} to {new_file_path}")
        else:
            print(f"No mapping found for {file_path}, Prefix: {prefix}")

# Rename both _1 and _2 files
rename_files('1')
rename_files('2')

print('Renaming process completed.')