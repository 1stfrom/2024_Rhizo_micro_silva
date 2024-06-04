#differ HN/LN for asv_fil
# Loop through each column in the first row
asv_fil_match <- asv_fil
for (i in seq_along(asv_fil_match[1,])) {
  # Modify the first row entries to keep only the first four characters
  asv_fil_match[, 1] <- substr(asv_fil_match[, 1], start = 1, stop = 4)
}

# Load the data from the text file
data <- read.table("/Users/yma18/Documents/GitHub/2024_Rhizo_micro_silva/data/2022_Root_Phenotyping_Data_checked.txt", header = TRUE, sep = "\t")  # Adjust sep based on your file's delimiter

# Check the third column and create datasets based on the condition
id_list_HN <- data[data[,3] == "HN", 1]  # Assuming 'HN' is in the third column and IDs in the first
id_list_LN <- data[data[,3] == "LN", 1]  # Assuming 'LN' is in the third column and IDs in the first

asv_fil_HN <- asv_fil_match[asv_fil_match[, 1] %in% id_list_HN, ]
asv_fil_LN <- asv_fil_match[asv_fil_match[, 1] %in% id_list_LN, ]


ps.clean.HN <- phyloseq(otu_table(asv_fil_HN),
               #sample_data(metadata_df),
               tax_table(tax_fil))
ps.clean.HN