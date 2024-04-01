library(dada2)
# File parsing
filtpathF <- "/Users/nathanma/Documents/PHDLife/mike_proj/data/largedata/F/filtered"
filtpathR <- "/Users/nathanma/Documents/PHDLife/mike_proj/data/largedata/R/filtered"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) #samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
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
saveRDS(seqtab, "/Users/nathanma/Documents/PHDLife/mike_proj/out/seqtab.rds")

#test protocal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


#####################################################
#Ver2 with debug
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  
  # Debugging print to check the number of merged sequences
  cat(paste("Number of sequences after merge for", sam, ":", nrow(merger)), "\n")
  
  mergers[[sam]] <- merger
}

# Additional check to ensure all elements in 'mergers' have the correct structure
if(any(sapply(mergers, function(x) is.null(x) || nrow(x) == 0))) {
  stop("Some merger elements are NULL or empty.")
}

# Now you can proceed with makeSequenceTable
seqtab <- makeSequenceTable(mergers)


####ver3 with debug
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  
  derepF <- derepFastq(filtFs[[sam]])
  if (is.null(derepF)) {
    cat("dereplication of forward reads failed for sample:", sam, "\n")
    next
  }
  
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  if (is.null(ddF)) {
    cat("DADA on forward reads failed for sample:", sam, "\n")
    next
  }
  
  derepR <- derepFastq(filtRs[[sam]])
  if (is.null(derepR)) {
    cat("dereplication of reverse reads failed for sample:", sam, "\n")
    next
  }
  
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  if (is.null(ddR)) {
    cat("DADA on reverse reads failed for sample:", sam, "\n")
    next
  }
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  if (is.null(merger) || nrow(merger) == 0) {
    cat("No sequences were merged for sample:", sam, "\n")
    next
  }
  
  cat("Merged sequences for sample:", sam, ", total:", nrow(merger), "\n")
  mergers[[sam]] <- merger
}

# After processing all samples, check the structure of the first few mergers
if (length(mergers) > 0) {
  cat("Checking the structure of the first merger element\n")
  str(mergers[[1]])
}

# Attempt to make a sequence table
cat("Attempting to create a sequence table...\n")
seqtab <- makeSequenceTable(mergers)

# Check if seqtab has been successfully created
if (!is.null(seqtab)) {
  cat("Sequence table created successfully. Dimensions:", dim(seqtab), "\n")
} else {
  cat("Failed to create sequence table. Please check the structure of mergers list.\n")
}
