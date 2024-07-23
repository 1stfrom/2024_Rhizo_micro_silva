library(vegan)
library(Ropt)
library(data.table)
library(phyloseq)
library(DESeq2)
library(tidyverse)
library(ape)
library("picante")
library("rstatix")

seqtab.nochim=readRDS("seqtab.nochim.rds")
taxa=readRDS("taxa.rds")
colnames(seqtab.nochim) <- paste0("ASV_", seq(1:ncol(seqtab.nochim)))
rownames(taxa) <- paste0("ASV_", seq(1:nrow(taxa)))
row.names(seqtab.nochim) <- substr(row.names(seqtab.nochim), 1, 4)
sampledata = read.csv("#phenotype_data/##.csv", header = FALSE)
colnames(sampledata) <- sampledata[1, ]
sampledata = sampledata[-1,]
rownames(sampledata) <- sampledata[ ,1]
sampledata = sampledata[ ,-1]

#Combine ASV,taxa into phyloseq 
OTU = otu_table(seqtab.nochim, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               tax_table(taxa))
#Phylo_tree
TRE = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

ps = merge_phyloseq(ps, sample_data(sampledata), TRE)
#flitration
ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("mitochondria"))
ps.clean
saveRDS(ps.clean,file="ps.clean.rds")

ps.clean.p0 <- filter_taxa(ps.clean, function(x) {sum(x > 3) >= (0.1*length(x))}, prune=TRUE)
ps.clean.p0

asv_fil = otu_table(ps.clean.p0)
asv2=asv_fil[apply(asv_fil, 1, sum)>1000,]
asv22=asv2+1
asv3=asv22@.Data
###ASV normalization
asv3=varianceStabilizingTransformation(asv3,blind = F)
ps.clean.re <- transform_sample_counts(ps.clean.p0, function(x) x / sum(x))
#ps.clean.re2 <- filter_taxa(ps.clean.re, function(x) sd(x) /mean(x)>3,TRUE)
cv=apply(asv3, 2, function(x)sd(x) /mean(x))
x=which(cv>1)
asv3=asv3[,x]
asv_fil=data.frame(ID=rownames(asv3),asv3)
asv_fil[,1]=unlist(strsplit(asv_fil[,1],split="_"))[seq(1,nrow(asv_fil)*2,by=2)]
asv_fil=asv_fil[asv_fil[,1]!="MC",]
tax_fil=tax_table(ps.clean.p0)
tax_fil=data.frame(ID=rownames(tax_fil),tax_fil)
write.table(asv_fil,"ASV_filter_DEseq_normalized.txt",sep="\t",quote=F,col.names = T,row.names = F)
write.table(tax_fil,"Taxa_filter.txt",sep="\t",quote=F,col.names = T,row.names = F)

