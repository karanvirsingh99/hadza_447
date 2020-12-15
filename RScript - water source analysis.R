#MAKING PCoA PLOTS - 26.11.20

# Loading packages
library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(ape)

# Helper functions
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
# Importing all data files
biom_file <- import_biom("./exported/table-with-taxonomy.biom")
metadata <- import_qiime_sample_data("tanzania_metadata.txt")
tree <- read_tree_greengenes("./exported/tree.nwk")
tree <- multi2di(tree)

# Combine all objects into phyloseq object
metadata$age_years = as.numeric(metadata$age_years)
physeq <- merge_phyloseq(biom_file, metadata, tree)

# Set set of random numbers
set.seed(711)

# taxonomic rank names
rank_names(physeq)
# Rename column names for taxonomic ranks
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus",
                                 "Species")
rank_names(physeq)

# Filter data based on metadata category with ==
human_samples <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                  sample_type == "feces" & age_years >= 3)

# Beta diversity PCoA plot
# Diversity requires rarefied taxa tables
# We know from Qiime that we should rarefy to 1250 so right palm
# samples level out
physeq_rar <- rarefy_even_depth(human_samples, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis with ordinate() and setting the
# method = PCoA, setting distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis and so on

#ordinatemultiple <- function(x) {
ord <- ordinate(physeq_rar_RA, method = "PCoA",
                distance = "Unweighted Unifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord,
                type = "sample",
                color = "water_source",
                title = "PCoA (Unweighted Unifrac) - no < 3yo") +
  # Adding text to plot
  #annotate(geom = "text",
          # label = "My great label",
          # x = - 0.08,
           #y = 0.18,
           #size = 4) + 
  # Manually adjust colours for points
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green", "brown"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED", "NOT COLLECTED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("unwUnifrac_water_no_below3yo.pdf", device = pdf, width = 30, height = 20, units = "cm")


##ORIGINAL GRAPH BASED ON WATER SOURCE, NO FILTERING
human_samples_total <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                  sample_type == "feces")



physeq_rar <- rarefy_even_depth(human_samples_total, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

ord_tot <- ordinate(physeq_rar_RA, method = "PCoA",
                distance = "Unweighted Unifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_tot,
                type = "sample",
                color = "water_source",
                title = "PCoA (Unweighted Unifrac) - all samples") +
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green", "brown"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED", "NOT COLLECTED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("unwUnifrac_water.pdf", device = pdf, width = 30, height = 20, units = "cm")

#FILTER OUT "NOT COLLECTED" WATER SOURCES
human_samples_filter_notcollected <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                  sample_type == "feces" & age_years >= 3 & water_source != "not_collected")

physeq_rar <- rarefy_even_depth(human_samples_filter_notcollected, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_filtered <- ordinate(physeq_rar_RA, method = "PCoA",
                distance = "Unweighted Unifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_filtered,
                type = "sample",
                color = "water_source",
                title = "PCoA (Unweighted Unifrac) - excluded < 3yo and not_collected") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("unwUnifrac_water_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Sex differentiated
human_samples_sex <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                                    )

physeq_rar <- rarefy_even_depth(human_samples_sex, sample.size = 14250)
  

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_sex <- ordinate(physeq_rar_RA, method = "PCoA",
                         distance = "Unweighted Unifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_sex,
                type = "sample",
                color = "water_source",
                shape = "sex",
                title = "PCoA (Unweighted Unifrac) - filtered by sex") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(labels=c("Female", "Male"),
                     values = c(1, 4)) +
  #stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source"),
         shape = guide_legend("sex")) + 
  theme_bw(base_size = 14)

ggsave("unwUnifrac_water_sex_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")
  
#Males ONLY
human_samples_males <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                      & sex != "FEMALE")

physeq_rar <- rarefy_even_depth(human_samples_males, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_males <- ordinate(physeq_rar_RA, method = "PCoA",
                    distance = "Unweighted Unifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_males,
                type = "sample",
                color = "water_source",
                title = "PCoA (Unweighted Unifrac) - males only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("unwUnifrac_water_males_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Females ONLY:
human_samples_females <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                        sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                      & sex != "MALE")

physeq_rar <- rarefy_even_depth(human_samples_females, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_females <- ordinate(physeq_rar_RA, method = "PCoA",
                      distance = "Unweighted Unifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_females,
                type = "sample",
                color = "water_source",
                title = "PCoA (Unweighted Unifrac) - females only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("unwUnifrac_water_females_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Season filtering:
human_samples <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                          sample_type == "feces" & age_years >= 3 &
                                          water_source != "not_collected")
                                        

physeq_rar <- rarefy_even_depth(human_samples, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_season <- ordinate(physeq_rar_RA, method = "PCoA",
                        distance = "Unweighted Unifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_season,
                type = "sample",
                color = "water_source",
                shape = "season",
                title = "PCoA (unweighted Unifrac) - seasons") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  #stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source"), shape = guide_legend("season")) + 
  theme_bw(base_size = 14)

ggsave("unwUnifrac_water_season_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")
