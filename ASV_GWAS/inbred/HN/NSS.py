import pandas as pd

# Step 1: Load the data files
data_blue = pd.read_csv('data_blue_HN_inb_836.txt', sep="\t")
nsslist = pd.read_csv('../../../SSlist.txt', sep="\t")

# Step 2: Extract the 'id' column (assumed to be the first column) from NSSlist
matching_ids = nsslist['id'].tolist()

# Step 3: Filter the rows in data_blue where the 'id' matches with those in NSSlist
matched_data = data_blue[data_blue['id'].isin(matching_ids)]

# Step 4: Save the matched data if needed
matched_data.to_csv('inb_HN_SS.txt', sep='\t', index=False)