library(dada2)
# Merge multiple runs (if necessary)
st1 <- readRDS("/Users/nathanma/Documents/PHDLife/mike_proj/out/seqtab.rds")
#st2 <- readRDS("path/to/run2/output/seqtab.rds")
#st3 <- readRDS("path/to/run3/output/seqtab.rds")
#st.all <- mergeSequenceTables(st1, st2, st3)

# Remove chimeras
seqtab <- removeBimeraDenovo(st1, method="consensus", multithread=TRUE)
#seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/Users/nathanma/Documents/PHDLife/mike_proj/database/silva_nr99_v138.1_train_set.fa", multithread=TRUE)

taxa<- assignTaxonomy(seqtab.nochim, "/Users/nathanma/Documents/PHDLife/mike_proj/database/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/nathanma/Documents/PHDLife/mike_proj/database/silva_species_assignment_v138.1.fa")
# Write to disk
saveRDS(seqtab, "/Users/nathanma/Documents/PHDLife/mike_proj/out/seqtab_final.rds") #sequence table saved
saveRDS(tax, "/Users/nathanma/Documents/PHDLife/mike_proj/out/tax_final.rds")