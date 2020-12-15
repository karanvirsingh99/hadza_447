##### ANIMAL ANALYSIS #####

library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(ape)
library(ggrepel)
library(ggpubr)


# Importing all data files
biom_file <- import_biom("table-with-taxonomy.biom")
metadata <- import_qiime_sample_data("tanzania_metadata.txt")
tree <- read_tree_greengenes("tree.nwk")
tree <- multi2di(tree)

# Combine all objects into phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

# Set set of random numbers
set.seed(711)

# Rename column names for taxonomic ranks
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus",
                                 "Species")
rank_names(physeq)

# Change "age_years" column to numeric

sample_data(physeq)$age_years <- as.numeric(sample_data(physeq)$age_years)
sample_data(physeq)$collection_timestamp <- as.Date(sample_data(physeq)$collection_timestamp,
                                                    format="%Y-%m-%d")

### ALL HUMAN AND ANIMAL SAMPLES (447_KS_7) ###

# Filter data based on metadata category
#all animal and human samples
physeq_ANIMALandHUMAN <- subset_samples(physeq, 
                                        env_feature %in% c("animal-associated habitat",
                                                           "human-associated habitat") &
                                          sample_type == "feces")


ANIMALandHUMAN_rar <- rarefy_even_depth(physeq_ANIMALandHUMAN, sample.size = 10435)


ANIMALandHUMAN_rar_RA <- transform_sample_counts(ANIMALandHUMAN_rar, function(x) x/sum(x))


ordinatemultiple <- function(p) {
  ord <- list()
  ord[["wunifrac"]] <- ordinate(p, method = "PCoA",
                                distance = "wunifrac")
  ord[["unifrac"]] <- ordinate(p, method = "PCoA",
                               distance = "unifrac")
  ord[["bray"]] <- ordinate(p, method = "PCoA",
                            distance = "bray")
  ord[["jaccard"]] <- ordinate(p, method = "PCoA",
                               distance = "jaccard")
  return(ord)
}


ord.animalANDhuman <- ordinatemultiple(ANIMALandHUMAN_rar_RA)

pdf("ALL PCOAS, ANIMAL+HUMAN.pdf", width=7, height=5)
for (x in c("wunifrac", "unifrac", "bray", "jaccard")) {
  plot <- plot_ordination(ANIMALandHUMAN_rar_RA,
                          ord.animalANDhuman[[x]],
                          type = "samples",
                          color = "env_feature",
                          title = paste0("PCoA ", x))+
    geom_point(size=3)+
    guides(color = guide_legend("SAMPLE TYPE"))+
    scale_color_discrete(labels=c("ANIMAL", "HUMAN"))
  print(plot)
}

pdf("ALL PCOAS, ANIMAL+HUMAN, with animal tags.pdf", width=7, height=5)
for (x in c("wunifrac", "unifrac", "bray", "jaccard")) {
  plot <- plot_ordination(ANIMALandHUMAN_rar_RA,
                          ord.animalANDhuman[[x]],
                          type = "samples",
                          color = "env_feature",
                          title = paste0("PCoA ", x))+
    geom_point(size=3)+
    guides(color = guide_legend("SAMPLE TYPE"))+
    scale_color_discrete(labels=c("ANIMAL", "HUMAN"))+
    geom_label_repel(aes(label = ifelse(sample_source =="VERVET_MONKEY",sample_source,"")),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')
  print(plot)
}

### ONLY 2 BUSH CAMPS, >=3YRS OLD (447_KS_8) ###

# Metadata based filtering in two steps (human and animal), then merge

#only human, hukamako and sengeli, 3yrs and older
physeq_HUMAN_twobush_3yrs <- subset_samples(physeq, 
                                            env_feature == "human-associated habitat" &
                                              bush_camp %in% c("HUKAMAKO", "SENGELI") &
                                              sample_type == "feces" &
                                              age_years >= 3)
#only animal samples
physeq_ANIMAL <- subset_samples(physeq, 
                                env_feature == "animal-associated habitat" &
                                  sample_type == "feces")

#human and animal, hukamako and sengeli, 3yrs and older
physeq_HUMAN_ANIMAL_2BUSH_3YRS <- merge_phyloseq(physeq_HUMAN_twobush_3yrs,
                                                 physeq_ANIMAL)

# Beta diversity PCoA plot
# Diversity requires rarefied taxa tables

HUMAN_ANIMAL_2BUSH_3YRS_rar <- rarefy_even_depth(physeq_HUMAN_ANIMAL_2BUSH_3YRS, sample.size = 10435)

# Convert to RA (relative abundance)
HUMAN_ANIMAL_2BUSH_3YRS_rar_RA <- transform_sample_counts(HUMAN_ANIMAL_2BUSH_3YRS_rar, function(x) x/sum(x))

ord.humanANDanimal2bush3yrs <- ordinatemultiple(HUMAN_ANIMAL_2BUSH_3YRS_rar_RA)


# Plot PCoA, with all non-human samples labelled

# Breaks for the PCoA timestamp
breaks.v2 <- seq(as.Date("2013/05/01"), as.Date("2014/10/1"), c("2 months"))

pdf("ALL PCOA HUMAN+ANIMAL, 3YRS2BUSH, with timestamp.pdf", useDingbats = FALSE,
    width=10, height=7)
for (x in c("wunifrac", "unifrac", "bray", "jaccard")) {
  
  plot <- plot_ordination(HUMAN_ANIMAL_2BUSH_3YRS_rar_RA,
                          ord.humanANDanimal2bush3yrs[[x]],
                          type = "samples",
                          color = "collection_timestamp",
                          shape= "env_feature")+
    # ggtitle(paste("PCoA", "unifrac"), subtitle="Human samples (HUKAMAKO and SENGELI, 3yrs and above) + animal samples")+
    geom_point(size=3)+
    scale_colour_gradientn(colours=c("#FFFFBF", 
                                     "#FDAE61",
                                     "#D53E4F",
                                     "#ABDDA4",
                                     "#66C2A5", 
                                     "#3288BD",
                                     "#FFFFBF", 
                                     "#FDAE61",
                                     "#D53E4F", 
                                     "#D53E4F"),
                           breaks=breaks.v2)+
    labs(colour="COLLECTION DATE", 
         shape="TYPE OF SAMPLE")+
    scale_shape(labels=c("ANIMAL", "HUMAN"))+
    geom_label_repel(aes(label = ifelse(!(sample_source =="HUMAN"),sample_source,"")),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')
  
  print(plot)
}
### BAR GRAPHS - THREE HUMANS AND TWO ANIMALS (447_KS_9) ###

#three outlier human samples in 2bush3yrs object
specialsamples<- c("11358.1589.fastq.gz", 	
                   "11358.1615.fastq.gz", 	
                   "11358.1592.fastq.gz")

#the two animal samples that cluser with the three outliers above
twoanimals <- c("11358.1609.fastq.gz", "11358.1173.fastq.gz")

# consistent colour scheme
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

twoanimalsandthreehumans <- HUMAN_ANIMAL_2BUSH_3YRS_rar %>%
  subset_samples(X.SampleID %in% c(twoanimals, specialsamples)) %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  #Transform to rel. abundance
  tax_glom(taxrank = "Phylum") %>%                      # agglomerate at phylum level
  psmelt() %>%                                          # Melt to long format
  filter(Abundance > 0.01) %>%                      #Filter low abundance OTUs
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

ggplot(twoanimalsandthreehumans, aes(x = Sample, y = Abundance, fill = Phylum)) +
  facet_grid(~env_feature, scale="free_x")+
  geom_bar(stat = "identity")+
  scale_x_discrete(labels=c("IMPALA", "HYRAX", "HUMAN", "HUMAN", "HUMAN"))+
  scale_fill_manual(values = twentyeightcolours)

### BAR GRAPHS - OUTLIER SAMPLES VS OTHER HUMAN SAMPLES (447_KS_10) ###

# Select the outlier samples
special <- subset_samples(HUMAN_2BUSH_3YRS_rar, X.SampleID %in% specialsamples)

# Take out all of the outlier samples
HUMAN_2BUSH_3YRS_rar_RA_NOSPECIAL <- HUMAN_2BUSH_3YRS_rar_RA %>% 
  subset_samples(!(X.SampleID %in% specialsamples))

SPECIALbar <- special %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  #Transform to rel. abundance
  tax_glom(taxrank = "Phylum") %>%                      # agglomerate at phylum level
  psmelt() %>%                                          # Melt to long format
  filter(Abundance > 0.01) %>%                      #Filter low abundance OTUs
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

NOTSPECIALbar.2014EW <- HUMAN_2BUSH_3YRS_rar_RA_NOSPECIAL %>%
  subset_samples(season == "2014-EW") %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Phylum)

NOTSPECIALbar<- HUMAN_2BUSH_3YRS_rar_RA_NOSPECIAL %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Phylum)

plot1 <- ggplot(SPECIALbar, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = twentyeightcolours)+
  scale_x_discrete(labels=NULL)

plot2 <- ggplot(NOTSPECIALbar.2014EW, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity")+
  scale_x_discrete(labels=NULL)+
  scale_fill_manual(values = twentyeightcolours)

plot3 <- ggplot(NOTSPECIALbar, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity")+
  scale_x_discrete(labels=NULL)+
  scale_fill_manual(values = twentyeightcolours)

ggarrange(
  plot1, plot2, plot3, labels = c("A", "B", "C"),
  common.legend = FALSE, legend = "right"
)


### BAR GRAPHS - MONKEYS AND HUMANS (447_KS_11) ###

vervet_monkeys_bar <- HUMAN_ANIMAL_2BUSH_3YRS_rar %>%
  subset_samples(sample_source %in% c("VERVET_MONKEY", "HUMAN")) %>% #only monkeys and humans
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  #Transform to rel. abundance
  tax_glom(taxrank = "Phylum") %>%                      # agglomerate at phylum level
  psmelt() %>%                                          # Melt to long format
  filter(Abundance > 0.01) %>%                      #Filter low abundance OTUs
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

plot <- ggplot(vervet_monkeys_bar, aes(x = Sample, y = Abundance, fill = Phylum)) +
  facet_grid(~host_common_name, scale="free_x")+
  geom_bar(stat = "identity")+
  labs(fill="PHYLUM")+
  scale_x_discrete(labels=NULL)+
  scale_fill_manual(values = twentyeightcolours)
# ggtitle("Taxa bar plots for vervet monkey samples and human samples",
#         subtitle="Only human samples from HUKAMAKO and SENGELI, Abundance >1%")
