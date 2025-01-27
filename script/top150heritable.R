HNi=read.table("../Heritability_HN_inb_836.txt",header = T,sep="\t")
LNi=read.table("../Heritability_LN_inb_836.txt",header = T,sep="\t")

HNi150 <- HNi[order(HNi$Heritability, decreasing = TRUE), ][1:150, ]
LNi150 <- LNi[order(LNi$Heritability, decreasing = TRUE), ][1:150, ]

overlap_traits <- intersect(HNi150$Trait, LNi150$Trait)
print(overlap_traits)
num_overlap <- length(overlap_traits)
print(paste("Number of overlapping traits:", num_overlap))

HNi_overlap <- HNi150[HNi150$Trait %in% overlap_traits, ]
LNi_overlap <- LNi150[LNi150$Trait %in% overlap_traits, ]

