library("tidyverse")
library("dada2")
library("gridExtra")
library("data.table")
no_of_cores=4
# File parsing
pathF <- "/path_to_forward_seq/F"
pathR <- "/path_to_reverse_seq/R"
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
filtpathF <- "/path_to_forward_seq/F/filtered"
filtpathR <- "/path_to_reverse_seq/R/filtered"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
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
saveRDS(seqtab, "/##/seqtab.rds") #seq table

####Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#View(seqtab.nochim)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim,file="seqtab.nochim.rds")
asv=as.data.frame(seqtab.nochim)
asv=data.frame(ID=rownames(asv),asv)
fwrite(asv,"/##/ASV.txt",sep="\t",quote=F,col.names = T,row.names = F)
