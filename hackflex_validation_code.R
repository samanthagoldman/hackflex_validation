---
title: "code for Hackflex library preparation enables low-cost metagenomic profiling manuscript"
output: html_document
date: "2023-01-18"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

```


```{r libraries}
library(vegan)
library(ggplot2)
library(phyloseq)
library(ape)
library(qiime2R)
library(tidyverse)
library(dplyr)
library(phytools)
library(scales)
library(cowplot)
library(ggforce)
library(reshape)
library(ggpubr)

```

```{r Figure1}

## Figure 1A ##

#feature table and metadata downloaded from qiita and imported into QIIME2 where it was filtered to only Zymo mock DNA community samples and then loaded into R

SVs_zymo<-read_qza("woltka_species_feature-table_zymoonly.qza")$data


taxasums_zymo<-summarize_taxa(SVs_zymo, taxonomy)$Genus


zymo_barplot<- taxa_barplot(taxasums_zymo, metadata, "dilution") +
  labs(y="") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

zymo_barplot 

#expected zymo barplot created manually using known abundances

## Figure 1B ##
bracken <- read.table("/Users/samgoldman/Documents/hackflex_validation/bracken_hackflex_validation.txt", header = TRUE, sep = "\t")


# Multiply all values by 100
bracken[, 2:ncol(bracken)] <- bracken[, 2:ncol(bracken)] * 100


bracken_plot <- pivot_longer(bracken, cols = c(2:19), names_to="sample_names", values_to = "abundance")

bracken_plot 


library(ggplot2)

species_abund_zymo_brack <- ggplot(bracken_plot, aes(x = reorder(Species, abundance), y = abundance)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.15)) +
  geom_hline(yintercept = 12, linetype = 2) +
  geom_hline(yintercept = 2, linetype = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels by 45 degrees

species_abund_zymo_brack <- ggplot(bracken_plot, aes(x = reorder(Species, -abundance), y = abundance)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.15)) +
  geom_hline(yintercept = 12, linetype = 2) +
  geom_hline(yintercept = 2, linetype = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.y = element_line(linewidth = 0.5, color = "black"),
        axis.text.y = element_text(margin = margin(r = 10)),
        panel.grid.major.y = element_line(linewidth = 0.5, color = "gray90"),
        panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(bracken_plot$abundance), by = 5))

 species_abund_zymo_brack


print(species_abund_zymo_brack)



species_abund_zymo_brack

# Create the stacked barplot

library(ggplot2)

# ... (your previous code)

# Create the boxplot with confidence intervals using stat_summary
species_abund_zymo_brack <- ggplot(bracken_plot, aes(x = reorder(Species, -abundance), y = abundance)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.15)) +
  geom_hline(yintercept = 12, linetype = 2) +
  geom_hline(yintercept = 2, linetype = 2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_line(linewidth = 0.5, color = "black"),
    axis.text.y = element_text(margin = margin(r = 10)),
    panel.grid.major.y = element_line(linewidth = 0.5, color = "gray90"),
    panel.grid.minor.y = element_blank()
  ) +
  scale_y_continuous(breaks = seq(0, max(bracken_plot$abundance), by = 5)) +
  stat_summary(
    fun.data = function(x) {
      y <- median(x)
      ymin <- quantile(x, 0.025)
      ymax <- quantile(x, 0.975)
      data.frame(y = y, ymin = ymin, ymax = ymax)
    },
    geom = "errorbar",
    width = 0.5,
    position = position_dodge(width = 0.75)
  )

species_abund_zymo_brack



# Melt the data frame to long format
bracken_long <- tidyr::pivot_longer(bracken, cols = -Species, names_to = "Sample", values_to = "Abundance")


bracken_barplot <- ggplot(bracken_long, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Abundance", fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bracken_barplot

write.table(bracken_long, "/Users/samgoldman/Documents/hackflex_validation/bracken_long_dilution.txt", sep = "\t")

bracken_long <- read.table("/Users/samgoldman/Documents/hackflex_validation/bracken_long_dilution.txt", header=TRUE)


# Convert Species column to factor with the custom order
bracken_long$Species <- factor(bracken_long$Species, levels = taxa_order)



# Define the custom order of species
taxa_order <- c(
  "Limosilactobacillus_fermentum", "Staphylococcus_aureus", "Listeria_monocytogenes",
  "Enterococcus_faecalis", "Bacillus_subtilis", "Escherichia_coli", "Salmonella_enterica",
  "Pseudomonas_aeruginosa", "Cryptococcus_neoformans", "Saccharomyces_cerevisiae"
)

# Assuming you have already read the data into bracken_long
bracken_long$Species <- factor(bracken_long$Species, levels = taxa_order)

# Reorder the data frame based on Species and Dilution
bracken_long <- bracken_long %>%
  arrange(Species, Dilution)

# Create the bar plot
bracken_barplot <- ggplot(bracken_long, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Abundance", fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~Dilution, scale = "free")

bracken_barplot




## Figure 1C ##

#MIQ score table created from running MIQ score program 

miq_scatter <- ggplot(MIQ_table, aes(x = Dilution, y = Score, color = Dilution)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

lm_model <- lm(Score ~ Dilution, data = MIQ_table)
print(summary(lm_model))

miq_scatter

```


```{r Figure2}

## Figure 2A ##

#feature table and metadata downloaded from qiita and imported into QIIME2 where it was filtered to only mouse samples and then loaded into R

metadata<-read_q2metadata("hackflexprep_allsamples.txt")

SVs<-read_qza("woltka_species_feature-table_miceonly.qza")$data

taxonomy<-read_qza("woltka_species_taxonomy.qza")$data %>% parse_taxonomy()

taxasums<-summarize_taxa(SVs, taxonomy)$Genus

mouse_barplot<- taxa_barplot(taxasums, metadata, "mousevmouse") #mousevmouse is a metadata column 
mouse_barplot

##Figure 2B ##

#DEICODE implementation of Robust Aitchison Distance calculated in QIIME2

mice_deicode_matrix <- read.table("distance-matrix_mice_deicode.tsv", header = TRUE)

m <- as.matrix(mice_deicode_matrix)
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2) <- c("sample1", "sample2", "distance")
 
# Join m2 and metadata based on sample1
m2 <- merge(m2, metadata[c("SampleID", "host_subject_id", "sequencing_meth")], by.x = "sample1", by.y = "SampleID")

# Rename columns for sample1
names(m2)[4:5] <- c("mouseid_sample1", "library_sample1")

m2$sample2 <- gsub("X", "", m2$sample2)


# Join m2 and metadata based on sample2
m2 <- merge(m2, metadata[c("SampleID", "host_subject_id", "sequencing_meth")], by.x = "sample2", by.y = "SampleID")
names(m2)[6:7] <- c("mouseid_sample2", "library_sample2")

#add group for comparison of all samples

m2$group <- ifelse(m2$mouseid_sample1 == m2$mouseid_sample2 & 
                   m2$library_sample1 == m2$library_sample2, 
                   "same_mouse_same_lib", 
                   ifelse(m2$mouseid_sample1 == m2$mouseid_sample2 & 
                          m2$library_sample1 != m2$library_sample2, 
                          "same_mouse_diff_lib", 
                          ifelse(m2$library_sample1 == "hackflex" & 
                                 m2$library_sample2 == "hackflex" & 
                                 m2$mouseid_sample1 != m2$mouseid_sample2, 
                                 "diff_mouse_hackflex_lib",
                                 ifelse(m2$library_sample1 == "tru" & 
                                        m2$library_sample2 == "tru" & 
                                        m2$mouseid_sample1 != m2$mouseid_sample2, 
                                        "diff_mouse_tru_lib",
                                        "all_diff"))))



m2_filtered <- m2 %>%
  filter(group != "all_diff")

group_order <- c("same_mouse_same_lib", "same_mouse_diff_lib", "diff_mouse_hackflex_lib", "diff_mouse_tru_lib")

m2_filtered$group <- factor(m2_filtered$group, levels = group_order)


#with jitter
p <- ggplot(m2_filtered, aes(x = group, y = distance)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_x_discrete(labels = c("Same Mouse Same Lib" = "same_mouse_same_lib",
                              "Same Mouse Diff Lib" = "same_mouse_diff_lib",
                              "Diff Mouse Hackflex Lib" = "diff_mouse_hackflex_lib",
                              "Diff Mouse Tru Lib" = "diff_mouse_tru_lib")) +
  labs(x = "Group", y = "Distance", title = "Distance between Samples by Group") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

p <- p + stat_compare_means(comparisons = list(c("same_mouse_same_lib", "same_mouse_diff_lib"),
                                               c("same_mouse_same_lib", "diff_mouse_hackflex_lib"),
                                               c("same_mouse_same_lib", "diff_mouse_tru_lib"),
                                               c("same_mouse_diff_lib", "diff_mouse_hackflex_lib"),
                                               c("same_mouse_diff_lib", "diff_mouse_tru_lib"),
                                               c("diff_mouse_hackflex_lib", "diff_mouse_tru_lib")),
                            label = "p.format", 
                            method = "wilcox.test",
                            hide.ns = TRUE)



# Subset data to include only relevant columns
data_for_test <- subset(m2_filtered, select = c("group", "distance"))

# Pairwise Wilcoxon rank sum test and extract p-values
pairwise_p_values <- pairwise.wilcox.test(data_for_test$distance, data_for_test$group, p.adjust.method = "bonferroni")$p.value

# Convert p-values to a matrix
p_value_matrix <- matrix(pairwise_p_values, nrow = length(levels(data_for_test$group)), ncol = length(levels(data_for_test$group)))
rownames(p_value_matrix) <- colnames(p_value_matrix) <- levels(data_for_test$group)

# Print p-value matrix and plot
p_value_matrix

p

## Figure 2C ##

#DEICODE implementation of Robust Aitchison Distance calculated in QIIME2

deicode <- read_qza("hackflexprep_mice_species_ordination_deicode.qza")


pcoa_mice <- deicode$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) 

pcoa_mice <- as_tibble(pcoa_mice)

pcoa_mice <- ggplot(pcoa_mice, aes(x=PC1, y=PC2, color=`host_subject_id`, group=`host_subject_id`, shape=`sample_zymoonly`)) +
  geom_mark_ellipse(aes(fill=`host_subject_id`), alpha=0.2) +
  geom_point(alpha=0.5) +
  theme_q2r() +
  geom_point(size=3) +
  scale_color_discrete(name="host_subject_id") +
  scale_fill_discrete(name="host_subject_id")

pcoa_mice

# Get percentage variance explained by each PC
variance <- deicode$data$ProportionExplained
variance <- paste0(round(variance * 100, 2), "%")

# Create PCoA plot with percentage explained by each PC
pcoa_mice <- deicode$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = host_subject_id, group = host_subject_id, shape = sample_zymoonly)) +
  geom_mark_ellipse() +
  geom_point(alpha = 0.5) +
  theme_q2r() +
  geom_point(size = 3) +
  scale_color_discrete(name = "host_subject_id") +
  labs(title = "PCoA Plot",
       x = paste0("PC1 (", variance[1], ")"),
       y = paste0("PC2 (", variance[2], ")"))


pcoa_mice


### Supplementary Figure 2 ###

#bracken_edit_GC.csv is the long form of the bracken output that has the GC values added as a column for each sample
bracken_gc <- read.table(bracken_edit_GC.csv, header = TRUE, sep = ",")

# Create scatter plot
gc_line_bracken <- ggplot(bracken_gc, aes(x = GC, y = delta_abund_fromknown)) +
  geom_point() +
  geom_smooth(method = lm)


print(gc_line_bracken)

# Fit a linear regression model
linear_model <- lm(delta_abund_fromknown ~ GC, data = bracken_gc)

summary_linear_model <- summary(linear_model)


print(summary_linear_model)


```