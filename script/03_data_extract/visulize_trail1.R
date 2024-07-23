#visulize
library(ggplot2) # Ensure ggplot2 is loaded for theme adjustment

png("Phylum_Relative abundance.png", height = 200, width = 400, units = "mm", res = 600)
par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
# Generate the plot with phyloseq::plot_bar and adjust legend text size
plot <- phyloseq::plot_bar(ps.clean.re, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") + 
  theme(legend.text = element_text(size = 8)) +
  theme(axis.title.x = element_text(size = 1))
print(plot)
dev.off()


##########################
library(phyloseq)
library(ggplot2)

# Begin PNG device
png("Phylum_Relative_Abundance.png", height = 300, width = 600, units = "mm", res = 1200)

# Generate the plot with phyloseq::plot_bar and adjust
plot <- phyloseq::plot_bar(ps.clean.re, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") + 
  theme(legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 0.5)) # Optional: Adjust x-axis text

# Print the plot
print(plot)

# Turn off the device
dev.off()

###########heat map
total = median(sample_sums(ps.clean.re))
ps.clean.re_abund <- filter_taxa(ps.clean.p0, function(x) sum(x > total*0.20) > 0, TRUE)


png("ASV_heat_map.png", height = 300, width = 600, units = "mm", res = 1200)
plot.heat <- phyloseq::plot_heatmap(ps.clean.re, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
                                    taxa.label = "Phylum", taxa.order = "Phylum", 
                                    trans=NULL, low="beige", high="red", na.value="beige")
print(plot)
dev.off()

#############tree
library("ape")

ps.clean.re_filtered <- prune_samples(sample_sums(ps) > 0, ps.clean.re)

random_tree = rtree(ntaxa(ps.clean.re_filtered), rooted=TRUE, tip.label=taxa_names(ps.clean.re_filtered))
plot(random_tree)
ps.clean.ret1 = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa), random_tree)
ps.clean.ret1

plot_tree(ps.clean.ret1, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(ps.clean.ret1, color="Depth", shape="Location", label.tips="taxa_names", ladderize="right", plot.margin=0.3)

#write.table(ps.clean.re_filtered,"tree.txt",sep="\t",quote=F,col.names = T,row.names = F)


ps.clean.re_filtered0 <- filter_taxa(ps.clean.re_filtered, function(x) sum(x > total*0.20) > 0, TRUE)
png("plot_heatmap.png", height = 300, width = 600, units = "mm", res = 1200)
plot.heat <- phyloseq::plot_heatmap(ps.clean.re_filtered, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
                                    taxa.label = "Phylum", taxa.order = "Phylum", 
                                    trans=NULL, low="beige", high="red", na.value="beige")
print(plot.heat)
dev.off()


plot_heatmap(ps.clean.re_filtered, taxa.label="Phylum")
