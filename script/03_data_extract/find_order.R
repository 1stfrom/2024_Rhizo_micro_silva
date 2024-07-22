# foucus on class Gammaproteobacteria
ps.class <- subset_taxa(ps, Class == "Gammaproteobacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!is.na(Order)) %>%
  subset_taxa(!is.na(Class)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("mitochondria"))
ps.class
saveRDS(ps.class,file="ps.class.rds")

ps.class.p0 <- filter_taxa(ps.class, function(x) {sum(x > 3) >= (0.1*length(x))}, prune=TRUE)
ps.class.p0 #taxa
asv_fil.class = otu_table(ps.class.p0)
asv2.class=asv_fil[apply(asv_fil.class, 1, sum)>1000,]
asv22.class=asv2.class+1
asv3.class=asv22.class@.Data
###ASV normalization
asv3.class=varianceStabilizingTransformation(asv3.class,blind = F)
ps.class.re <- transform_sample_counts(ps.class.p0, function(x) x / sum(x))
#ps.class.re2 <- filter_taxa(ps.class.re, function(x) sd(x) /mean(x)>3,TRUE)
cv.class=apply(asv3.class, 2, function(x)sd(x) /mean(x))
x=which(cv.class>1)
asv3.class=asv3.class[,x]
asv_fil.class=data.frame(ID=rownames(asv3.class),asv3.class)
asv_fil.class[,1]=unlist(strsplit(asv_fil.class[,1],split="_"))[seq(1,nrow(asv_fil.class)*2,by=2)]
asv_fil.class=asv_fil.class[asv_fil.class[,1]!="MC",]
tax_fil.class=tax_table(ps.class.p0)
tax_fil.class=data.frame(ID=rownames(tax_fil.class),tax_fil.class)
write.table(asv_fil.class,"ASV_pro_normalized.txt",sep="\t",quote=F,col.names = T,row.names = F)
write.table(tax_fil.class,"Taxa_pro_filter.txt",sep="\t",quote=F,col.names = T,row.names = F)

png("Class_Relative abundance.png",height = 100,width =210,units = "mm",res=600)
par(mfrow=c(1,1),mar=c(4,4,1,1))
phyloseq::plot_bar(ps.class.re, fill = "Order") +
  geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack")
dev.off()