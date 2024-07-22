import os

# Specify the directory containing the files you want to rename
#directory_path = '/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/2023_root_micro/F/filtered'

# Iterate over each file in the directory
for filename in os.listdir(directory_path):
    # Check if the filename contains an underscore
    if "_" in filename:
        # Construct the new filename by replacing underscores with nothing
        new_filename = filename.replace("_", "")
        # Construct the full old and new file paths
        old_file_path = os.path.join(directory_path, filename)
        new_file_path = os.path.join(directory_path, new_filename)
        # Rename the file
        os.rename(old_file_path, new_file_path)
        print(f'Renamed "{filename}" to "{new_filename}"')

print("All files have been renamed.")

import os
# Set the directory containing the files
directory_path = '/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/2023_root_micro/F/filtered'

# Iterate over each file in the directory
for filename in os.listdir(directory_path):
    # Check if the file ends with "R2001.fastq.gz"
    if filename.endswith("R2001.fastq.gz") and not filename.endswith("_R2001.fastq.gz"):
        # Construct the new filename by inserting an underscore before "R2001"
        new_filename = filename.replace("R2001.fastq.gz", "_R2001.fastq.gz")
        # Construct the full old and new file paths
        old_file_path = os.path.join(directory_path, filename)
        new_file_path = os.path.join(directory_path, new_filename)
        # Rename the file
        os.rename(old_file_path, new_file_path)
        print(f'Renamed "{filename}" to "{new_filename}"')

print("File renaming completed.")
