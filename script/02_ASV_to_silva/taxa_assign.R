library("tidyverse")
library("dada2")
library("gridExtra")
library("phyloseq")
library("Biostrings")
library("data.table")
library("DESeq2")
###Taxonomy Inference
seqtab.nochim=readRDS("seqtab.nochim.rds")
taxa=readRDS("taxa.rds")
gc(full=TRUE)
silva <- "/path_to_database/silva_nr99_v138.1_train_set.fa"
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=TRUE, verbose=TRUE)
saveRDS(taxa,file="taxa.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
seqtab.taxa.plus=cbind('#seq'=rownames(taxa),t(seqtab.nochim),taxa)
write.table(seqtab.taxa.plus,"ASV.taxon.species.txt",sep="\t",quote=F,col.names = T,row.names = F)