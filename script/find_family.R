#foucus on family Proteobacteria
ps.pro <- subset_taxa(ps, Phylum == "Proteobacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!is.na(Class)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("mitochondria"))
ps.pro
saveRDS(ps.pro,file="ps.pro.rds")

ps.pro.p0 <- filter_taxa(ps.pro, function(x) {sum(x > 3) >= (0.1*length(x))}, prune=TRUE)
ps.pro.p0 #taxa
asv_fil.pro = otu_table(ps.pro.p0)
asv2.pro=asv_fil[apply(asv_fil.pro, 1, sum)>1000,]
asv22.pro=asv2.pro+1
asv3.pro=asv22.pro@.Data
###ASV normalization
asv3.pro=varianceStabilizingTransformation(asv3.pro,blind = F)
ps.pro.re <- transform_sample_counts(ps.pro.p0, function(x) x / sum(x))
#ps.pro.re2 <- filter_taxa(ps.pro.re, function(x) sd(x) /mean(x)>3,TRUE)
cv.pro=apply(asv3.pro, 2, function(x)sd(x) /mean(x))
x=which(cv.pro>1)
asv3.pro=asv3.pro[,x]
asv_fil.pro=data.frame(ID=rownames(asv3.pro),asv3.pro)
asv_fil.pro[,1]=unlist(strsplit(asv_fil.pro[,1],split="_"))[seq(1,nrow(asv_fil.pro)*2,by=2)]
asv_fil.pro=asv_fil.pro[asv_fil.pro[,1]!="MC",]
tax_fil.pro=tax_table(ps.pro.p0)
tax_fil.pro=data.frame(ID=rownames(tax_fil.pro),tax_fil.pro)
write.table(asv_fil.pro,"ASV_pro_normalized.txt",sep="\t",quote=F,col.names = T,row.names = F)
write.table(tax_fil.pro,"Taxa_pro_filter.txt",sep="\t",quote=F,col.names = T,row.names = F)

png("Pro_Relative abundance.png",height = 100,width =210,units = "mm",res=600)
par(mfrow=c(1,1),mar=c(4,4,1,1))
phyloseq::plot_bar(ps.pro.re, fill = "Class") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack")
dev.off()










