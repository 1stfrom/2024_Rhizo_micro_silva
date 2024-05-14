library(vegan)
library(Ropt)
library(data.table)
library(phyloseq)
library(DESeq2)
library(tidyverse)
library(ape)

library("picante")
library("rstatix")

seqtab.nochim=readRDS("out/seqtab.nochim.rds")
taxa=readRDS("out/taxa.rds")
colnames(seqtab.nochim) <- paste0("ASV_", seq(1:ncol(seqtab.nochim)))
rownames(taxa) <- paste0("ASV_", seq(1:nrow(taxa)))
row.names(seqtab.nochim) <- substr(row.names(seqtab.nochim), 1, 4)
sampledata = read.csv("data/2022_Root_Phenotyping_Data_checked.csv", header = FALSE)
colnames(sampledata) <- sampledata[1, ]
sampledata = sampledata[-1,]
rownames(sampledata) <- sampledata[ ,1]
sampledata = sampledata[ ,-1]


OTU = otu_table(seqtab.nochim, taxa_are_rows = TRUE)
TAX = tax_table(taxa)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               tax_table(taxa))
TRE = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

ps = merge_phyloseq(ps, sample_data(sampledata), TRE)

ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("mitochondria"))
ps.clean

ps.clean.pg <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 25}, prune=TRUE)
ps.clean.pg
ps.clean.p8 <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= (0.1*length(x))}, prune=TRUE)
ps.clean.p8
ps.clean.p0 <- filter_taxa(ps.clean, function(x) {sum(x > 3) >= (0.1*length(x))}, prune=TRUE)
ps.clean.p0

ps.clean.p125 <- filter_taxa(ps.clean, function(x) {sum(x > 10) >= (0.125*length(x))}, prune=TRUE)
ps.clean.p125
ps.clean.p75 <- filter_taxa(ps.clean, function(x) {sum(x > 3) >= (0.075*length(x))}, prune=TRUE)
ps.clean.p75
ps.clean.p5 <- filter_taxa(ps.clean, function(x) {sum(x > 3) >= (0.05*length(x))}, prune=TRUE)
ps.clean.p5
ps.clean.p25 <- filter_taxa(ps.clean, function(x) {sum(x > 3) >= (0.025*length(x))}, prune=TRUE)
ps.clean.p25




ps.clean.p3 <- filter_taxa(ps.clean, function(x) {sum(x > 7) >= (0.005*length(x))}, prune=TRUE)
ps.clean.p4 <- filter_taxa(ps.clean, function(x) {sum(x > 7) >= (0.015*length(x))}, prune=TRUE)
ps.clean.p3
ps.clean.p4


ps.clean.re <- transform_sample_counts(ps.clean.p0, function(x) x / sum(x))
FSfr = filter_taxa(ps.clean.re, function(x) sum(x) > .005, TRUE)
ps.clean.re125 <- transform_sample_counts(ps.clean.p125, function(x) x / sum(x))


#Absolute abundance (stacked bar plots)
png("Phylum_Relative_Abundance_HNLN.png", height = 300, width = 600, units = "mm", res = 1200)

phyloseq::plot_bar(ps.clean.p0, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  facet_wrap(~trt, scales = "free_x")
  theme(legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 10), 
      axis.text.x = element_text(angle = 45, hjust = 1, size = 0.5))

dev.off()



# Agglomerate to phylum-level and rename
ps.clean.phylum = phyloseq::tax_glom(ps.clean.p0, taxrank = rank_names(ps.clean.p0)[2])
phyloseq::taxa_names(ps.clean.phylum) = phyloseq::tax_table(ps.clean.phylum)[,"Phylum"]
otu_table(ps.clean.phylum)[1:5,1:5]


#Melt and plot
png("Phylum_Relative_Abundance_each_indi.png", height = 300, width = 600, units = "mm", res = 1200)
phyloseq::psmelt(ps.clean.phylum) %>%
  ggplot(data = ., aes(x = trt, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free", nrow = 2) +
  ggtitle("Phylum absolute abundance by trt")

dev.off()

# Agglomerate to phylum-level and rename
png("Phylum_Relative_Abundance_each_geno.png", height = 300, width = 600, units = "mm", res = 1200)
phyloseq::psmelt(ps.clean.phylum) %>%
  ggplot(data = ., aes(x = geno, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free", nrow = 2) +
  ggtitle("Phylum absolute abundance by trt") +
  theme(axis.text.x = element_text(size=0.2))

dev.off()




##Relative abundance
ps.clean.re <- transform_sample_counts(ps.clean.p0, function(x) x / sum(x))
png("Phylum_Relative_Abundance.png", height = 300, width = 600, units = "mm", res = 1200)
phyloseq::plot_bar(ps.clean.re, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  facet_wrap(~trt, scales = "free_x")
dev.off()



# alpha diversity should be calculated before filtering on abundance and prevalence
tree = phyloseq::phy_tree(ps)
samp = data.frame(phyloseq::otu_table(ps))

#Generate a data.frame with adiv measures
adiv <- data.frame(
  phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps)), geno, trt)) %>%
  dplyr::select(-se.chao1)

measures = colnames(adiv)[1:7]
#head(adiv)




####among geno diversity
adiv %>%
  group_by( geno ) %>%
  dplyr::summarise(median_Observed = median(Observed),
                   median_Chao1 = median(Chao1),
                   median_Shannon = median(Shannon),
                   median_Simpson = median(Simpson),
                   median_InvSimpson = median(InvSimpson),
                   median_Fisher = median(Fisher),
                   median_PD = median(PD))
####among trt diversity
adiv %>%
  group_by( trt ) %>%
  dplyr::summarise(median_Observed = median(Observed),
                   median_Chao1 = median(Chao1),
                   median_Shannon = median(Shannon),
                   median_Simpson = median(Simpson),
                   median_InvSimpson = median(InvSimpson),
                   median_Fisher = median(Fisher),
                   median_PD = median(PD))


#Plot adiv measures
adiv %>%
  gather(key = metric, value = value, measures) %>% #c("Observed", "Shannon", "Chao1", "Simpson", "PD")
  mutate(metric = factor(metric, levels = measures)) %>% #c("Observed", "Shannon", "Chao1", "Simpson", "PD")
  ggplot(aes(x = geno, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = geno), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measures") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") 


res = adiv %>%
  gather(key = metric, value = value, measures) %>% 
  mutate(metric = factor(metric, levels = measures)) %>%
  group_by(metric) %>%
  rstatix::wilcox_test( value ~ trt, p.adjust.method = "none") %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "trt", scales = "free_y")






ps.beta = ps.clean.re

sample_sums <- rowSums(otu_table(ps.beta))
print(names(sample_sums[sample_sums == 0]))
ps.beta <- prune_samples(sample_sums > 0, ps.beta)

# 检查OTU表中的无穷大或NA值
otu <- otu_table(ps.beta)
if (any(is.infinite(otu)) || any(is.na(otu))) {
  # 将无穷大值或NA替换为0（根据情况可以选择其他替换值）
  otu[is.infinite(otu) | is.na(otu)] <- 0
  otu_table(ps.beta) <- otu
}




ordBC <- ordinate(ps.beta, "PCoA", "bray")
#str(ordBC)
ordJC <- ordinate(ps.beta, "PCoA", "jaccard")
ordUF <- ordinate(ps.beta, "PCoA", "unifrac")
ordwUF <- ordinate(ps.beta, "PCoA", "wunifrac")



smpID <- sample_data(ps.beta)$geno

# Keep first 2 vectors (latent variables, PCs) of each distance matrix
df <- rbind(data.frame(ordBC$vectors[,1:2], geno = smpID,method = 'BC'),
            data.frame(ordJC$vectors[,1:2], geno = smpID,method = 'Jaccard'),
            #data.frame(ordUF$vectors[,1:2], geno = smpID,method = 'unifrac'),
            data.frame(ordwUF$vectors[,1:2], geno = smpID,method = 'wunifrac'))

# add sample_data info
df <- merge(df, data.frame(sample_data(ps.beta)), by = 'geno')




ggplot(data = df, aes(Axis.1,Axis.2, 
                      color = trt) ) + 
  geom_point() + 
  stat_ellipse() + 
  facet_wrap(~method,scales = 'free')











asv_fil = otu_table(ps.clean.p9)
asv2=asv_fil[apply(asv_fil, 1, sum)>1000,]
asv22=asv2+1
asv3=asv22@.Data
###ASV normalization
asv3=varianceStabilizingTransformation(asv3,blind = F)
ps.clean.re <- transform_sample_counts(ps.clean.p9, function(x) x / sum(x))
#ps.clean.re2 <- filter_taxa(ps.clean.re, function(x) sd(x) /mean(x)>3,TRUE)
cv=apply(asv3, 2, function(x)sd(x) /mean(x))
x=which(cv>1)
asv3=asv3[,x]
asv_fil=data.frame(ID=rownames(asv3),asv3)
asv_fil[,1]=unlist(strsplit(asv_fil[,1],split="_"))[seq(1,nrow(asv_fil)*2,by=2)]
asv_fil=asv_fil[asv_fil[,1]!="MC",]
tax_fil=tax_table(ps.clean.p9)
tax_fil=data.frame(ID=rownames(tax_fil),tax_fil)
write.table(asv_fil,"ASV_filter_DEseq_normalized.txt",sep="\t",quote=F,col.names = T,row.names = F)
write.table(tax_fil,"Taxa_filter.txt",sep="\t",quote=F,col.names = T,row.names = F)





ps.re <- transform_sample_counts(ps.clean.p0, function(x) x / sum(x))

####get genus table
gen <- tax_glom(ps,taxrank='Genus')
gen1 <- subset_taxa(gen, Genus != "uncultured") 
taxa_names(gen1) <- as.data.frame(as(tax_table(gen1),'matrix'))$Genus # Rename taxa using Genus names
gen2 <- filter_taxa(gen1, function (x) {sum(x > 0) >= 832}, prune=TRUE)
genus=t(as.data.frame(otu_table(gen2)))
genus=data.frame(name=rownames(genus),genus)

cv=apply(genus[,-1], 2, function(x)sd(x) /mean(x))
hist(cv,breaks = 50,xlab="",main="CV of ASV across the pop")
x=which(cv>2)
genus=genus[,c(1,x+1)]
genus2=decostand(genus[,-1], method="normalize",MARGIN=1)
genus3=cbind(genus[,1],genus2)
colnames(genus3)[1]="ID"
write.table(genus3, file = "Normalized_Genus_high_abundance3_CV2.txt", sep = "\t", quote = F,row.names = F,col.name=T)

# Calculate relative abundance of each genus:
ge_re <- data.frame(taxa_sums(gen1)/sum(taxa_sums(gen1)))


#####################Diversity
# alpha diversity should be calculated before filtering on abundance and prevalence
# plot rarefaction curve
otu=otu_table(ps)

class(otu) <- "matrix" # as.matrix() will do nothing
otu <- t(otu)
rarecurve(otu, step=200, label = F, sample = c(1000,2000))

rarecurve(t(otu_table(ps)), step=200, label = F, sample = c(1000,2000))
# rarefaction 
ps.rare = rarefy_even_depth(ps, 1000, rngseed = 777777, replace = FALSE)
# generate a data.frame with alpha diversity measures
adiv = data.frame(
  "Observed" = phyloseq::estimate_richness(ps.rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps.rare, measures = "Shannon"))
adiv = cbind(adiv, sample_data(ps.rare))
write.table(adiv, file = "alpha_diversity_table_1000.txt", sep = "\t", quote = F)


# beta diversity should be calculated using relative abundance
ps.clean.p0 <- filter_taxa(ps, function (x) {sum(x > 0) >= 106.1}, prune=TRUE)








ps.clean.phylum = phyloseq::tax_glom(ps.clean.p0, taxrank = rank_names(ps.clean.p0)[2])
phyloseq::taxa_names(ps.clean.phylum) = phyloseq::tax_table(ps.clean.phylum)[,"Phylum"]
otu_table(ps.clean.phylum)[1:5,1:5]





phyloseq::psmelt(ps.clean.phylum) %>%
  ggplot(data = ., aes(y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free", nrow = 2) +
  ggtitle("Phylum absolute abundance by time")




