import pandas as pd

# Step 1: Load the data files
inb_HN_SS = pd.read_csv('inb_HN_NSS.txt', sep="\t")
geno_list = pd.read_csv('geno_list.txt', sep="\t")

# Step 2: Extract the 'id' column from geno_list and inb_HN_SS
geno_ids = geno_list['id'].tolist()

# Step 3: Merge the data
# Perform a left join, keeping all ids from geno_list and assigning NaN for those not in inb_HN_SS
ready_for_GWAS = pd.merge(geno_list, inb_HN_SS, on='id', how='left')

# Step 4: Save the result to a new file as ready_for_GWAS.txt
ready_for_GWAS.to_csv('inb_HN_NSS_GWAS.txt', sep='\t', index=False)