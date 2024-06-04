library("tidyverse")
library("dada2")
library("gridExtra")
library("data.table")
no_of_cores=4
# File parsing
pathF <- "/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/F"
pathR <- "/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/R"
filtpathF <- file.path(pathF, "filtered")
filtpathR <- file.path(pathR, "filtered")
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(280,220), rm.phix=TRUE,
              trimLeft=c(20, 20), maxEE=c(2,2),
              compress=TRUE, verbose=TRUE, multithread=no_of_cores)

# File parsing
filtpathF <- "/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/F/filtered"
filtpathR <- "/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/R/filtered"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
#sample.names <- gsub("_clean_R1.fastq.gz|/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/F/filtered/","",filtFs)
#sample.namesR <- gsub("_clean_R2.fastq.gz|/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/data/raw_reads/R/filtered/","",filtRs)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, verbose=TRUE, multithread=no_of_cores)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, verbose=TRUE, multithread=no_of_cores)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/out/seqtab.rds") #seq table

####Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#View(seqtab.nochim)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim,file="seqtab.nochim.rds")
asv=as.data.frame(seqtab.nochim)
asv=data.frame(ID=rownames(asv),asv)
fwrite(asv,"/Users/nathanma/Documents/PHDLife/mike_proj/out/ASV.txt",sep="\t",quote=F,col.names = T,row.names = F)
#seqtab=readRDS("seqtab.rds")

#######3
library("tidyverse")
library("dada2")
library("gridExtra")
library("phyloseq")
library("Biostrings")
library("data.table")
library("DESeq2")
###Taxonomy Inference
##seqtab.nochim=readRDS("seqtab.nochim.rds")
##taxa=readRDS("taxa.rds")
gc(full=TRUE)
silva <- "/Users/nathanma/Documents/PHDLife/mike_proj/database/silva_nr99_v138.1_train_set.fa"
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=TRUE, verbose=TRUE)
saveRDS(taxa,file="taxa.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
seqtab.taxa.plus=cbind('#seq'=rownames(taxa),t(seqtab.nochim),taxa)
write.table(seqtab.taxa.plus,"ASV.taxon.species.txt",sep="\t",quote=F,col.names = T,row.names = F)





###Phylogeny Reconstruction
##Before we can calculate a phylogenetic tree we need a FASTA file of the ASV sequences.
# Name the sequence file
#seqtab.nochim=readRDS("seqtab.nochim.rds")
sequence_outfile <- "dada2_nochim.fa"?
print(paste0('Writing FASTA file ', sequence_outfile, ' with ', ncol(seqtab.nochim), ' sequences.'))

file.remove(sequence_outfile) # Make sure the file does not exist
for (i in 1:ncol(seqtab.nochim)){
  cat(paste0(">ASV_", i), file = sequence_outfile, sep="\n", append = TRUE)
  cat(colnames(seqtab.nochim)[i], file = sequence_outfile, sep="\n", append= TRUE)
}

# We need to make the external packages MAFFT and FastTree available for our use;
# this is what the next three lines do.
source(file.path(Sys.getenv("LMOD_PKG"), "init/R"))
module("load", "mafft")
module("load", "fasttree")
#/work/jyanglab/nathanma/config/R/packages/mafft-linux64
# Multiple sequence alignment with MAFFT
system("mafft --auto --thread -1 dada2_nochim.fa > dada2_nochim_mafft_msa.fa", intern=TRUE)

# Phylogenetic tree reconstruction with FastTree
system("fasttree -gtr -nt < dada2_nochim_mafft_msa.fa > dada2_nochim.tree", intern=TRUE)

list.files(pattern = "*dada2*")






seqtab.nochim=readRDS("/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/out/2023_root_micro/seqtab.nochim.rds")
taxa=readRDS("/work/jyanglab/nathanma/project/16SRNA_vali/2024_vali_silva/out/2023_root_micro/taxa.rds")
seqtab.nochim.0 <- seqtab.nochim
taxa.0 <- taxa
##Change colnames(seqtab.nochim) and rownames(taxa) to the same: ASV_1, ASV_2…
colnames(seqtab.nochim) <- paste0("ASV_", seq(1:ncol(seqtab.nochim)))
rownames(taxa) <- paste0("ASV_", seq(1:nrow(taxa)))
###use phyloseq() to make the data object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
          #     sample_data(metadata_df),
               tax_table(taxa))
ps
###Step 1, remove Archaea, unknowns, chloroplasts;
ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("mitochondria"))
ps.clean
saveRDS(ps.clean,file="ps.clean.rds")
#ps.clean=readRDS("ps.clean.rds")
##Step 2, filter taxa (ASV) with low abundance (< 2);
ps.clean.p0 <- filter_taxa(ps.clean, function(x) {sum(x > 3) >= (0.1*length(x))}, prune=TRUE)
ps.clean.p0 #taxa
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

png("Phylum_Relative abundance.png",height = 100,width =210,units = "mm",res=600)
par(mfrow=c(1,1),mar=c(4,4,1,1))
phyloseq::plot_bar(ps.clean.re, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")
dev.off()