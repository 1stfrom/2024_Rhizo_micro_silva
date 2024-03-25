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

taxa <- assignTaxonomy(seqtab.nochim, "/Users/nathanma/Documents/PHDLife/mike_proj/database/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/nathanma/Documents/PHDLife/mike_proj/database/silva_species_assignment_v138.1.fa")
# Write to disk
saveRDS(seqtab, "/Users/nathanma/Documents/PHDLife/mike_proj/out/seqtab_final.rds") #sequence table saved
saveRDS(tax, "/Users/nathanma/Documents/PHDLife/mike_proj/out/tax_final.rds")

write.csv(tax,file="/Users/nathanma/Documents/PHDLife/mike_proj/out/test_protcal/taxa.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")

####
saveRDS(seqtab.nochim,file="seqtab.nochim.rds")
asv=as.data.frame(seqtab.nochim)
asv=data.frame(ID=rownames(asv),asv)
fwrite(asv,"/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/out/ASV.txt",sep="\t",quote=F,col.names = T,row.names = F)
#seqtab=readRDS("seqtab.rds")


gc(full=TRUE) 
silva <- "/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/database_99/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=no_of_cores, verbose=TRUE)
saveRDS(taxa,file="taxa.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
seqtab.taxa.plus=cbind('#seq'=rownames(taxa),t(seqtab.nochim),taxa)
write.table(seqtab.taxa.plus,"ASV.taxon.species.txt",sep="\t",quote=F,col.names = T,row.names = F)