#### Bush camps location + water source analysis ####

##########Script 1##########
#Generate PCoA plotted on Bray-Curtis, unweighted UniFrac,
#weighted UniFrac for for final filtered dataset that includes 
#all bush camps. 

# Loading packages
library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(ape)

#Importing all data files previously exported from QIIME
biom_file_1 <- import_biom("table-with-taxonomy-1.biom")
metadata <- import_qiime_sample_data("tanzania_metadata.txt")
tree <- read_tree_greengenes("tree.nwk")

#Convert Qiime's multichotomous tree into a dichotomous one
tree <- multi2di(tree)

#Combine all objects into phyloseq object
physeq <- merge_phyloseq(biom_file_1, metadata, tree)

#Set set of random numbers
set.seed(711)

#Look at taxonomic rank names
rank_names(physeq)

#Rename column names for taxonomic ranks
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus",
                                 "Species")
#Create table rarified to 11000
physeq_rar <- rarefy_even_depth(physeq, sample.size = 11000)

#Convert to relative abundance
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#Bray-curtis 
ord <- ordinate(physeq_rar_RA, method = "PCoA",
                distance = "bray")
#Plot data
plot_ordination(physeq_rar_RA,
                ord,
                type = "sample",
                color = "bush_camp",
                title = "",) +
  geom_point(size = 3)+
  guides(colour = guide_legend("BUSH CAMP"))+
  theme_bw(base_size = 14)

#Save the plot
ggsave("PCoA_bray_age_allcamps.pdf", width=10, height=8)

#Unweighted UniFrac
ord1 <- ordinate(physeq_rar_RA, method = "PCoA",
                 distance = "unifrac")
# Plot data
plot_ordination(physeq_rar_RA,
                ord1,
                type = "sample",
                color = "bush_camp",
                title = "",) +
  geom_point(size = 3)+
  guides(colour = guide_legend("BUSH CAMP"))+
  theme_bw(base_size = 14)

#Save the plot
ggsave("PCoA_unifrac_age_allcamps.pdf", width=10, height=8)

#Weighted UniFrac
ord2 <- ordinate(physeq_rar_RA, method = "PCoA",
                 distance = "wunifrac")
# Plot data
plot_ordination(physeq_rar_RA,
                ord2,
                type = "sample",
                color = "bush_camp",
                title = "",) +
  geom_point(size = 3)+
  guides(colour = guide_legend("BUSH CAMP"))+
  theme_bw(base_size = 14)

#Save the plot
ggsave("PCoA_wunifrac_age_allcamps.pdf", width=10, height=8)

#############################################################
###########Script 2##########
# Loading packages
library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(ape)

# Importing all data files
biom_file_2 <- import_biom("table-with-taxonomy-2.biom")
metadata <- import_qiime_sample_data("tanzania_metadata.txt")
tree <- read_tree_greengenes("tree.nwk")

#Convert Qiime's multichotomous tree into a dichotomous one
tree <- multi2di(tree)

#Combine all objects into phyloseq object
physeq <- merge_phyloseq(biom_file_2, metadata, tree)


#Set set of random numbers
set.seed(711)

#Look at taxonomic rank names
rank_names(physeq)

#Rename column names for taxonomic ranks
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus",
                                 "Species")
#Create table rarified to 11000
physeq_rar <- rarefy_even_depth(physeq, sample.size = 11000)

#Convert to relative abundance
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))


#Bray-curtis 
ord3 <- ordinate(physeq_rar_RA, method = "PCoA",
                 distance = "bray")

# Plot data
plot_ordination(physeq_rar_RA,
                ord3,
                type = "sample",
                color = "bush_camp",
                shape = "water_source",
                title = "") +
  geom_point(size = 2) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("BUSH CAMP")) + 
  guides(shape = guide_legend("WATER SOURCE"))+
  theme_bw(base_size = 14)

#Save the plot
ggsave("PCoA_Bray_final.pdf", width = 10, height = 8)



#Unweighted-unifrac
ord4 <- ordinate(physeq_rar_RA, method = "PCoA",
                 distance = "unifrac")

#Plot data
plot_ordination(physeq_rar_RA,
                ord4,
                type = "sample",
                color = "bush_camp",
                shape = "water_source",
                title = "") +
  geom_point(size = 2)+
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("BUSH CAMP")) + 
  guides(shape = guide_legend("WATER SOURCE"))+
  theme_bw(base_size = 14)

#Save the plot
ggsave("PCoA_unweighted_final.pdf", width =10, height =8)


#Weighted-unifrac
ord5 <- ordinate(physeq_rar_RA, method = "PCoA",
                 distance = "wunifrac")
#Plot data
plot_ordination(physeq_rar_RA,
                ord5,
                type = "sample",
                color = "bush_camp",
                shape = "water_source",
                title = "") +
  geom_point(size = 2)+
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("BUSH CAMP")) + 
  guides(shape = guide_legend("WATER SOURCE"))+
  theme_bw(base_size = 14)

#Save the plot
ggsave("PCoA_weighted_final.pdf", width =10, height =8)

##############################################################
##########Script 3##########
#Identify the outliers

#Look for outliers in unweighted unifrac
ordinationvectors2 <- ord4[["vectors"]]
#Save the outliers as a list
special<- c("11358.1589.fastq.gz", "11358.1592.fastq.gz", "11358.1267.fastq.gz", "11358.1241.fastq.gz", "11358.1615.fastq.gz")

#Look for outliers in weighted unifrac
ordinationvectors3 <- ord5[["vectors"]]
#Outliers in both unweighted and weighted unifrac are the same

#Subset the samples
special <- subset_samples(physeq_rar_RA, X.SampleID %in% special)

#Make the taxa bar plots
SPECIALbar <- special %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt() %>%
  filter(Abundance>0.01) %>%
  arrange(Phylum)

#Make the colour scheme based on Karanvir's in order to maintain consistency
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

twentyeightcolours <- sample(col_vector, 28)

twentyeightcolours <- c(p__ = "#FC8D62", `p__[Thermi]` = "#FDCDAC", p__Acidobacteria = "#FFFF99", 
                        p__Actinobacteria = "#CAB2D6", p__Armatimonadetes = "#CCCCCC", 
                        p__Bacteroidetes = "#F4CAE4", p__Chlorobi = "#E5D8BD", p__Chloroflexi = "#A6CEE3", 
                        p__Crenarchaeota = "#B3B3B3", p__Cyanobacteria = "#DECBE4", p__Elusimicrobia = "#BF5B17", 
                        p__Euryarchaeota = "#7FC97F", p__FBP = "#66C2A5", p__Fibrobacteres = "#CBD5E8", 
                        p__Firmicutes = "#B3E2CD", p__Fusobacteria = "#D9D9D9", p__Gemmatimonadetes = "#B3DE69", 
                        p__Lentisphaerae = "#FBB4AE", p__Nitrospirae = "#BEBADA", p__Planctomycetes = "#E78AC3", 
                        p__Proteobacteria = "#BC80BD", p__Spirochaetes = "#1B9E77", p__SR1 = "#A6761D", 
                        p__Synergistetes = "#999999", p__Tenericutes = "#E6AB02", p__TM7 = "#666666", 
                        p__Verrucomicrobia = "#E31A1C", `p__WPS-2` = "#8DA0CB")

names(twentyeightcolours) = c("p__", "p__[Thermi]", "p__Acidobacteria", "p__Actinobacteria", 
                              "p__Armatimonadetes", "p__Bacteroidetes", "p__Chlorobi", "p__Chloroflexi", 
                              "p__Crenarchaeota", "p__Cyanobacteria", "p__Elusimicrobia", "p__Euryarchaeota", 
                              "p__FBP", "p__Fibrobacteres", "p__Firmicutes", "p__Fusobacteria", 
                              "p__Gemmatimonadetes", "p__Lentisphaerae", "p__Nitrospirae", 
                              "p__Planctomycetes", "p__Proteobacteria", "p__Spirochaetes", 
                              "p__SR1", "p__Synergistetes", "p__Tenericutes", "p__TM7", "p__Verrucomicrobia", 
                              "p__WPS-2")
#Plot the graph
ggplot(SPECIALbar, aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = twentyeightcolours)

#Save the plot
ggsave("Taxa_barplots_final.pdf", width =10, height =8)
