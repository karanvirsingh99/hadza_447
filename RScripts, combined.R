#### WATER SOURCE ANALYSIS ####
### BRAY CURTIS PLOTS ###
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
                distance = "bray")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord,
                type = "sample",
                color = "water_source",
                title = "PCoA (Bray Curtis) - no < 3yo") +
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

ggsave("bray_water_no_below3yo.pdf", device = pdf, width = 30, height = 20, units = "cm")


##ORIGINAL GRAPH BASED ON WATER SOURCE, NO FILTERING
human_samples_total <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                        sample_type == "feces")



physeq_rar <- rarefy_even_depth(human_samples_total, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

ord_tot <- ordinate(physeq_rar_RA, method = "PCoA",
                    distance = "bray")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_tot,
                type = "sample",
                color = "water_source",
                title = "PCoA (Bray Curtis) - all samples") +
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green", "brown"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED", "NOT COLLECTED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("bray_water.pdf", device = pdf, width = 30, height = 20, units = "cm")

#FILTER OUT "NOT COLLECTED" WATER SOURCES
human_samples_filter_notcollected <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected")

physeq_rar <- rarefy_even_depth(human_samples_filter_notcollected, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_filtered <- ordinate(physeq_rar_RA, method = "PCoA",
                         distance = "bray")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_filtered,
                type = "sample",
                color = "water_source",
                title = "PCoA (Bray Curtis) - excluded < 3yo and not_collected") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("bray_water_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Sex differentiated
human_samples_sex <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
)

physeq_rar <- rarefy_even_depth(human_samples_sex, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_sex <- ordinate(physeq_rar_RA, method = "PCoA",
                    distance = "bray")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_sex,
                type = "sample",
                color = "water_source",
                shape = "sex",
                title = "PCoA (Bray Curtis) - filtered by sex") +
  
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

ggsave("bray_water_sex_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Males ONLY
human_samples_males <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                        sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                      & sex != "FEMALE")

physeq_rar <- rarefy_even_depth(human_samples_males, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_males <- ordinate(physeq_rar_RA, method = "PCoA",
                      distance = "bray")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_males,
                type = "sample",
                color = "water_source",
                title = "PCoA (Bray Curtis) - males only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("bray_water_males_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Females ONLY:
human_samples_females <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                          sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                        & sex != "MALE")

physeq_rar <- rarefy_even_depth(human_samples_females, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_females <- ordinate(physeq_rar_RA, method = "PCoA",
                        distance = "bray")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_females,
                type = "sample",
                color = "water_source",
                title = "PCoA (Bray Curtis) - females only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("bray_water_females_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Season filtering:
human_samples <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                  sample_type == "feces" & age_years >= 3 &
                                  water_source != "not_collected")


physeq_rar <- rarefy_even_depth(human_samples, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_season <- ordinate(physeq_rar_RA, method = "PCoA",
                       distance = "bray")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_season,
                type = "sample",
                color = "water_source",
                shape = "season",
                title = "PCoA (Bray Curtis) - seasons") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  #stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source"), shape = guide_legend("season")) + 
  theme_bw(base_size = 14)

ggsave("bray_water_season_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

### Jaccard Plot ###

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
                distance = "jaccard")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord,
                type = "sample",
                color = "water_source",
                title = "PCoA (Jaccard) - no < 3yo") +
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

ggsave("jaccard_water_no_below3yo.pdf", device = pdf, width = 30, height = 20, units = "cm")


##ORIGINAL GRAPH BASED ON WATER SOURCE, NO FILTERING
human_samples_total <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                        sample_type == "feces")



physeq_rar <- rarefy_even_depth(human_samples_total, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

ord_tot <- ordinate(physeq_rar_RA, method = "PCoA",
                    distance = "jaccard")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_tot,
                type = "sample",
                color = "water_source",
                title = "PCoA (Jaccard) - all samples") +
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green", "brown"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED", "NOT COLLECTED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("jaccard_water.pdf", device = pdf, width = 30, height = 20, units = "cm")

#FILTER OUT "NOT COLLECTED" WATER SOURCES
human_samples_filter_notcollected <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected")

physeq_rar <- rarefy_even_depth(human_samples_filter_notcollected, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_filtered <- ordinate(physeq_rar_RA, method = "PCoA",
                         distance = "jaccard")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_filtered,
                type = "sample",
                color = "water_source",
                title = "PCoA (Jaccard) - excluded < 3yo and not_collected") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("jaccard_water_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Sex differentiated
human_samples_sex <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
)

physeq_rar <- rarefy_even_depth(human_samples_sex, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_sex <- ordinate(physeq_rar_RA, method = "PCoA",
                    distance = "jaccard")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_sex,
                type = "sample",
                color = "water_source",
                shape = "sex",
                title = "PCoA (Jaccard) - filtered by sex") +
  
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

ggsave("jaccard_water_sex_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Males ONLY
human_samples_males <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                        sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                      & sex != "FEMALE")

physeq_rar <- rarefy_even_depth(human_samples_males, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_males <- ordinate(physeq_rar_RA, method = "PCoA",
                      distance = "jaccard")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_males,
                type = "sample",
                color = "water_source",
                title = "PCoA (Jaccard) - males only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("jaccard_water_males_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Females ONLY:
human_samples_females <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                          sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                        & sex != "MALE")

physeq_rar <- rarefy_even_depth(human_samples_females, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_females <- ordinate(physeq_rar_RA, method = "PCoA",
                        distance = "jaccard")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_females,
                type = "sample",
                color = "water_source",
                title = "PCoA (Jaccard) - females only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("jaccard_water_females_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Season filtering:
human_samples <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                  sample_type == "feces" & age_years >= 3 &
                                  water_source != "not_collected")


physeq_rar <- rarefy_even_depth(human_samples, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_season <- ordinate(physeq_rar_RA, method = "PCoA",
                       distance = "jaccard")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_season,
                type = "sample",
                color = "water_source",
                shape = "season",
                title = "PCoA (Jaccard) - seasons") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  #stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source"), shape = guide_legend("season")) + 
  theme_bw(base_size = 14)

ggsave("jaccard_water_season_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

### Unweighted UniFrac Plots ###

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

### Weighted UniFrac Analysis ###

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
                distance = "wunifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord,
                type = "sample",
                color = "water_source",
                title = "PCoA (Weighted Unifrac) - no < 3yo") +
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

ggsave("wunifrac_water_no_below3yo.pdf", device = pdf, width = 30, height = 20, units = "cm")


##ORIGINAL GRAPH BASED ON WATER SOURCE, NO FILTERING
human_samples_total <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                        sample_type == "feces")



physeq_rar <- rarefy_even_depth(human_samples_total, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

ord_tot <- ordinate(physeq_rar_RA, method = "PCoA",
                    distance = "wunifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_tot,
                type = "sample",
                color = "water_source",
                title = "PCoA (Weighted Unifrac) - all samples") +
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green", "brown"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED", "NOT COLLECTED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("wunifrac_water.pdf", device = pdf, width = 30, height = 20, units = "cm")

#FILTER OUT "NOT COLLECTED" WATER SOURCES
human_samples_filter_notcollected <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected")

physeq_rar <- rarefy_even_depth(human_samples_filter_notcollected, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_filtered <- ordinate(physeq_rar_RA, method = "PCoA",
                         distance = "wunifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_filtered,
                type = "sample",
                color = "water_source",
                title = "PCoA (Weighted Unifrac) - excluded < 3yo and not_collected") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("wunifrac_water_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Sex differentiated
human_samples_sex <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                      sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
)

physeq_rar <- rarefy_even_depth(human_samples_sex, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_sex <- ordinate(physeq_rar_RA, method = "PCoA",
                    distance = "wunifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_sex,
                type = "sample",
                color = "water_source",
                shape = "sex",
                title = "PCoA (Weighted Unifrac) - filtered by sex") +
  
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

ggsave("wunifrac_water_sex_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Males ONLY
human_samples_males <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                        sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                      & sex != "FEMALE")

physeq_rar <- rarefy_even_depth(human_samples_males, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_males <- ordinate(physeq_rar_RA, method = "PCoA",
                      distance = "wunifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_males,
                type = "sample",
                color = "water_source",
                title = "PCoA (Weighted Unifrac) - males only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("wunifrac_water_males_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Females ONLY:
human_samples_females <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                          sample_type == "feces" & age_years >= 3 & water_source != "not_collected"
                                        & sex != "MALE")

physeq_rar <- rarefy_even_depth(human_samples_females, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_females <- ordinate(physeq_rar_RA, method = "PCoA",
                        distance = "wunifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_females,
                type = "sample",
                color = "water_source",
                title = "PCoA (Weighted Unifrac) - females only") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  scale_shape_manual(values = "cross") +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source")) + 
  theme_bw(base_size = 14)

ggsave("wunifrac_water_females_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")

#Season filtering:
human_samples <- subset_samples(physeq, env_feature == "human-associated habitat" &
                                  sample_type == "feces" & age_years >= 3 &
                                  water_source != "not_collected")


physeq_rar <- rarefy_even_depth(human_samples, sample.size = 14250)


# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

#ordinatemultiple <- function(x) {
ord_season <- ordinate(physeq_rar_RA, method = "PCoA",
                       distance = "wunifrac")
#}
# Plot data
plot_ordination(physeq_rar_RA,
                ord_season,
                type = "sample",
                color = "water_source",
                shape = "season",
                title = "PCoA (Weighted Unifrac) - seasons") +
  
  scale_colour_manual(values = c("blue", "red",
                                 "black", "yellow", "green"),
                      labels = c("SPRING", "STREAM",
                                 "WELL",
                                 "DRY RIVER BED", "RIVER BED")) +
  #stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("water source"), shape = guide_legend("season")) + 
  theme_bw(base_size = 14)

ggsave("wunifrac_water_season_no_below3yo_notcollected.pdf", device = pdf, width = 30, height = 20, units = "cm")



#### Bush camps location + water source analysis ####

##Script 1##
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


###Script 2###
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

###Script 3###
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
