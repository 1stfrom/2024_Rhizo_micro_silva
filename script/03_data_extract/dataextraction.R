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

ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("mitochondria"))
ps.clean