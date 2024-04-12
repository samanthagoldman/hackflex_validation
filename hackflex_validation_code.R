---
title: "Hackflex validation code"
output: html_document
date: "2024-03-01"
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r, echo = FALSE}

```

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
library(ggforce)
library(reshape)
library(ggpubr)
library(compositions)
library(RVAideMemoire)
```


```{r Figure 1}

##Fig 1A##

# Kraken table created from https://ccb.jhu.edu/software/kraken2/

zymo_kraken <- read.table("zymo_combined_Kraken2_report_mpa_format.txt", header = TRUE, sep = "\t")

zymo_kraken <- as.data.frame(zymo_kraken)

kraken_zymo_long <- pivot_longer(zymo_kraken, cols = -species, names_to = "Sample", values_to = "Value")

# Normalize values to 100%
sample_sums <- aggregate(Value ~ Sample, data = kraken_zymo_long, sum)

# Merge sums with the original dataframe
kraken_zymo_long <- merge(kraken_zymo_long, sample_sums, by = "Sample", suffixes = c("", "_sum"))

# Calculate the normalized values
kraken_zymo_long$Normalized_Value <- (kraken_zymo_long$Value / kraken_zymo_long$Value_sum) * 100

# Make the barplot
kraken_zymo_normalized <- ggplot(kraken_zymo_long, aes(x = Sample, y = Normalized_Value, fill = species)) +
  geom_bar(stat = "identity") +
  labs(title = "Normalized Taxa Stacked Barplot",
       x = "Sample",
       y = "Percentage",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  kraken_zymo_normalized 


# Created the expected zymo barplot created manually using known abundances


## Figure 1B ##

kraken_boxplot <- ggplot(kraken_zymo_long, aes(x = Sample, y = Value, fill = species)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Abundance", fill = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Figure 1C ##

# MIQ_table created from https://github.com/Zymo-Research/miqScoreShotgunPublic

miq_scatter <- ggplot(MIQ_table, aes(x = Dilution, y = Score, color = Dilution)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

lm_model <- lm(Score ~ Dilution, data = MIQ_table)
print(summa
miq_scatter

```


```{r Figure 2}

## Figure 2A ##

# Feature table and metadata downloaded from Qiita and imported into QIIME2 where it was filtered and then loaded into R

metadata_nextera<-read_q2metadata("metadata.txt")

SVs_nextera<-read_qza("woltka_species_feature-table.qza")$data

taxonomy_nextera<-read_qza("/woltka_species_taxonomy_table.qza")$data %>% parse_taxonomy()

taxasums_nextera<-summarize_taxa(SVs_nextera, taxonomy_nextera)$Genus

mouse_barplot_nextera<- taxa_barplot(taxasums_nextera, metadata_nextera, "test")

mouse_barplot_nextera


## Figure 2B ##

# DEICODE implementation of Robust Aitchison Distance calculated in QIIME2

mice_deicode_matrix_nextera <- read.table("distance-matrix.tsv", header = TRUE)

# Find the row and column index containing "195474.15506.blank"
blank_col <- which(colnames(mice_deicode_matrix_nextera) == "X195474.15506.blank")


# Remove the row and column containing "blank"
mice_deicode_matrix_nextera <- mice_deicode_matrix_nextera[-8, -blank_col_index]

m_nextera <- as.matrix(mice_deicode_matrix_nextera)
m2_nextera <- melt(m_nextera)[melt(upper.tri(m_nextera))$value,]
names(m2_nextera) <- c("sample1", "sample2", "distance")
 
# Join m2 and metadata based on sample1
m2_nextera <- merge(m2_nextera, metadata_nextera[c("SampleID", "lib_meth", "mouse_id", "lib_dilution")], by.x = "sample1", by.y = "SampleID")


# Rename columns for sample1
names(m2_nextera)[4:5] <- c("library_sample1", "mouseid_sample1")

m2_nextera$sample2 <- gsub("X", "", m2_nextera$sample2)


# Join m2 and metadata based on sample2
m2_nextera <- merge(m2_nextera, metadata_nextera[c("SampleID", "lib_meth", "mouse_id", "lib_dilution")], by.x = "sample2", by.y = "SampleID")
names(m2_nextera)[7:8] <- c("library_sample2", "mouseid_sample2")

#add group for comparison of all samples
m2_nextera$group <- ifelse(m2_nextera$mouseid_sample1 == m2_nextera$mouseid_sample2 & 
                   m2_nextera$library_sample1 == m2_nextera$library_sample2, 
                   "same_mouse_same_lib", 
                   ifelse(m2_nextera$mouseid_sample1 == m2_nextera$mouseid_sample2 & 
                          m2_nextera$library_sample1 != m2_nextera$library_sample2, 
                          "same_mouse_diff_lib", 
                          ifelse(m2_nextera$library_sample1 == "hackflex" & 
                                 m2_nextera$library_sample2 == "hackflex" & 
                                 m2_nextera$mouseid_sample1 != m2_nextera$mouseid_sample2, 
                                 "diff_mouse_hackflex_lib",
                                 ifelse(m2_nextera$library_sample1 == "nextera" & 
                                        m2_nextera$library_sample2 == "nextera" & 
                                        m2_nextera$mouseid_sample1 != m2_nextera$mouseid_sample2, 
                                        "diff_mouse_nextera_lib",
                                        "all_diff"))))


# Filter out comparison that's unnecessary for figure
m2_filtered_nextera <- m2_nextera %>%
  filter(group != "all_diff")

group_order <- c("diff_mouse_nextera_lib", "diff_mouse_hackflex_lib", "same_mouse_diff_lib", "same_mouse_same_lib")

m2_filtered_nextera$group <- factor(m2_filtered_nextera$group, levels = group_order)


# Plot
p_nextera <- ggplot(m2_filtered_nextera, aes(x = group, y = distance)) +
  geom_boxplot() +
  geom_jitter(width = 0.15) +
  scale_x_discrete(labels = c("Same Mouse Same Lib" = "same_mouse_same_lib",
                              "Same Mouse Diff Lib" = "same_mouse_diff_lib",
                              "Diff Mouse Hackflex Lib" = "diff_mouse_hackflex_lib",
                              "Diff Mouse Nextera Lib" = "diff_mouse_tru_lib")) +
  labs(x = "Group", y = "Distance", title = "Distance between Samples by Group") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))


# Pairwise t test
read.delim("illuminadnaprep_deicode_withmd.txt", sep=" ")->dnaprep
pairwise.perm.t.test(dnaprep$distance,dnaprep$group)




## Figure 2C ##

# DEICODE implementation of Robust Aitchison Distance calculated in QIIME2

deicode_nextera <- read_qza("hackflex_species_ordination_deicode.qza")


pcoa_mice_nextera <- deicode_nextera$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata_nextera) 

pcoa_mice_nextera <- as_tibble(pcoa_mice_nextera)

# Get percentage variance explained by each PC
variance_nextera <- deicode_nextera$data$ProportionExplained
variance_nextera <- paste0(round(variance_nextera * 100, 2), "%")

# Create PCoA plot with percentage explained by each PC
pcoa_mice_nextera <- deicode_nextera$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata_nextera_pcoa) %>%
  ggplot(aes(x=PC1, y=PC2, color=`mouse_id`, group=`mouse_id`, shape=`lib_meth`)) +
  geom_mark_ellipse() +
  geom_point(alpha = 0.5) +
  theme_q2r() +
  geom_point(size = 3) +
  scale_color_discrete(name = "mouse_id") +
  labs(title = "PCoA Plot",
       x = paste0("PC1 (", variance_nextera[1], ")"),
       y = paste0("PC2 (", variance_nextera[2], ")"))


pcoa_mice_nextera

#PERMANOVA

adonis2(formula = mice_deicode_matrix_nextera ~ mouse_id * lib_meth, data = metadata_nextera)


adonis2(formula = mice_deicode_matrix_nextera ~ lib_meth * mouse_id, data = metadata_nextera)


```





```{r Supplementary Figures}


## Supplementary Figure 1 created with MIQ Score https://github.com/Zymo-Research/miqScoreShotgunPublic ##

##Supplementary Figure 2##

# First transpose the output from kraken, then load into R
zymo_kraken <- read.table("zymo_kraken_mpa_format.txt", header = TRUE, sep = "\t")

zymo_sample <- zymo_kraken$sample

expected <- read.table("expected_abund.txt", header = TRUE)

zymo_kraken[, -1] <- sapply(zymo_kraken[, -1], as.numeric)
expected[, -1] <- sapply(expected[, -1], as.numeric)

combined_zymo_kraken <- rbind(zymo_kraken[, -1],expected[1, -1])

#CLR transform
clr_combined_zymo_kraken <- Hotelling:::clr(combined_zymo_kraken)

# Add sample column with "expected" for the last entry
clr_combined_zymo_kraken <- cbind(sample = c(zymo_kraken$sample, "expected"), clr_combined_zymo_kraken)

clr_combined_zymo_kraken_differences<-clr_combined_zymo_kraken[,-1]-rbind(clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1],clr_combined_zymo_kraken[19,-1])

clr_combined_zymo_kraken_differences[-19,]

combined_clr_bind_diff_final <- cbind(clr_combined_zymo_kraken[-19,1], clr_combined_zymo_kraken_differences[-19,])

combined_clr_bind_diff_final <- t(combined_clr_bind_diff_final)

combined_clr_bind_diff_final <- as.data.frame(combined_clr_bind_diff_final)

# Set the first row as header
colnames(combined_clr_bind_diff_final) <- combined_clr_bind_diff_final[1, ]

# Remove the first row (since it's now the header)
combined_clr_bind_diff_final <- combined_clr_bind_diff_final[-1, ]

combined_clr_bind_diff_final <- t(combined_clr_bind_diff_final)

combined_clr_bind_diff_final <- as.data.frame(combined_clr_bind_diff_final)

zymo_kraken_clr_long <- pivot_longer(combined_clr_bind_diff_final, cols = colnames(combined_clr_bind_diff_final), names_to = "species", values_to = "clr_abundance_diff")

zymo_md <- read.table("zymo_community_metadata.txt", header = TRUE)

# Merge with metadata
zymo_kraken_clr_long <- merge(zymo_kraken_clr_long, zymo_md, by.x = "species", by.y = "species", all.x = TRUE)

# Calculate average CLR_diff for each Species

zymo_kraken_clr_long$clr_abundance_diff <- as.numeric(zymo_kraken_clr_long$clr_abundance_diff)

avg_forplotting <- aggregate(clr_abundance_diff ~ species, data = zymo_kraken_clr_long, FUN = mean)

avg_forplotting <- merge(avg_forplotting, additional_data, by.x = "species", by.y = "Species")

# Create scatter plot
gc_line_kraken <- ggplot(zymo_kraken_clr_long, aes(x = GC, y = clr_abundance_diff)) + 
geom_point() + 
  scale_color_manual(values = species_palette, drop = FALSE)
gc_line_kraken

print(gc_line_kraken)


# Add relevant columns to avg_forplotting
additional_data <- aggregate(cbind(type, genome_size, GC, gram) ~ Species, data = kraken_gc, FUN = function(x) x[1])
avg_forplotting <- merge(avg_forplotting, additional_data, by.x = "zymo_species", by.y = "Species")

avg_forplotting$GC <- as.numeric(avg_forplotting$GC)


# Add on average values
gc_line_kraken <- gc_line_kraken +
  geom_point(data = avg_forplotting, aes(x = GC, y = clr_abundance_diff), color = "red", size = 3)

# Plot
gc_line_kraken <- gc_line_kraken +
  labs(x = "% GC content", 
       y = "Change in Observed Relative Abundance Compared 
       to Known Abundance, CLR Transformed") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))


print(gc_line_kraken)



# Fit linear regression model
lm_model_gc <- lm(clr_abundance_diff ~ GC, data = avg_forplotting)

summary(lm_model_gc)

# Add linear regression line to the plot
gc_line_kraken<-gc_line_kraken + geom_smooth(data = avg_forplotting, method = "lm", se = FALSE, color = "blue")
gc_line_kraken


# Bias statistics

# Filter out yeasts
avg_kraken_bacteria <- avg_forplotting[avg_forplotting$gram != "yeast", ]


# Gram status
anova_model_gram <- aov(CLR_diff ~ gram, data = avg_bracken_bacteria)

summary(anova_model_gram)

# Domain

anova_model_type <- aov(clr_abundance_diff ~ type, data = avg_forplotting)

summary(anova_model_type)


# Genome size

lm_model_size <- lm(clr_abundance_diff ~ genome_size, data = avg_forplotting)

summary(lm_model_size)



#Genome size, bacteria only

lm_model_size_bact <- lm(clr_abundance_diff ~ genome_size, data = avg_kraken_bacteria)

summary(lm_model_size_bact)

#FDR correction of p values
FDR_corrected <- p.adjust(c(0.311, 0.0733, 0.2932, 0.5366, 0.003617), method="bonferroni")



##Supplementary Figure 3##

mouse_bias_taxa <- read.table("mouse_bias_taxa.txt", header = TRUE, sep = "\t")

mouse_sample <- mouse_bias_taxa$clade_name

mouse_bias_taxa[, -1] <- sapply(mouse_bias_taxa[, -1], as.numeric)

mouse_bias_taxa_clr <- clr_lite(mouse_bias_taxa[, -1], samples_are="rows", method="const", replicates = 1000)

rownames(mouse_bias_taxa_clr) <- mouse_bias_taxa[, 1] 

write.table(mouse_bias_taxa_clr, "mouse_bias_taxa_clr.txt", quote = FALSE, sep = "\t")

#Subtract differences and then reupload

mouse_bias_taxa_clr_diff <- read.table ("mouse_bias_taxa_clr_diff.txt", header = TRUE, sep = "\t")

mouse_bias_taxa_clr_diff <- as.data.frame(mouse_bias_taxa_clr_diff)

mouse_bias_taxa_clr_diff_long <- pivot_longer(mouse_bias_taxa_clr_diff, cols = -X, names_to = "species", values_to = "clr_abundance_diff")

zymo_kraken_clr_long <- pivot_longer(clr_combined_zymo_kraken, cols = -sample, names_to = "species", values_to = "clr_abundance_diff")

names(mouse_bias_taxa_clr_diff_long)[1] <- "sample"
names(mouse_bias_taxa_clr_diff_long)[2] <- "clade_name"

mouse_md <- read.table("mouse_bias_metadata.txt", header = TRUE)

#merge with metadata
mouse_bias_taxa_clr_diff_long <- merge(mouse_bias_taxa_clr_diff_long, mouse_md, by.x = "clade_name", by.y = "clade_name", all.x = TRUE)


# Calculate average CLR_diff for each Species

mouse_bias_taxa_clr_diff_long$clr_abundance_diff <- as.numeric(mouse_bias_taxa_clr_diff_long$clr_abundance_diff)

avg_forplotting_mouse <- aggregate(clr_abundance_diff ~ clade_name, data = mouse_bias_taxa_clr_diff_long, FUN = mean)

avg_forplotting_mouse <- merge(avg_forplotting_mouse, mouse_md, by.x = "clade_name", by.y = "clade_name")

# Create scatter plot
gc_line_mouse <- ggplot(mouse_bias_taxa_clr_diff_long, aes(x = GC_content, y = clr_abundance_diff)) + 
geom_point() 
gc_line_mouse

avg_forplotting_mouse$GC_content <- as.numeric(avg_forplotting_mouse$GC_content)

# Add average values
gc_line_mouse <- gc_line_mouse +
  geom_point(data = avg_forplotting_mouse, aes(x = GC_content, y = clr_abundance_diff), color = "red", size = 3)

gc_line_mouse  <- gc_line_mouse  +
  labs(x = "% GC content", 
       y = "Change in Hackflex Relative Abundance Compared 
       to TruSeq Relative Abundance, CLR Transformed") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))

gc_line_mouse

# Linear regression
lm_model_gc_mouse <- lm(clr_abundance_diff ~ GC_content, data = avg_forplotting_mouse)

summary(lm_model_gc_mouse)

# Add linear regression line to the plot
gc_line_mouse<-gc_line_mouse + geom_smooth(data = avg_forplotting_mouse, method = "lm", se = FALSE, color = "blue")
gc_line_mouse


# Gram status statistics

# Fit the two-way ANOVA model
anova_model_gram <- aov(clr_abundance_diff ~ gram_status, data = avg_forplotting_mouse)

summary(anova_model_gram)

#FDR correction
FDR_corrected <- p.adjust(c( 0.0985, 0.524), method="bonferroni")
FDR_corrected



## Supplementary Figure 4A ##

#feature table and metadata downloaded from qiita and imported into QIIME2 where it was filtered and then loaded into R

metadata<-read_q2metadata("metadata.txt")

SVs<-read_qza("woltka_species_feature-table.qza")$data

taxonomy<-read_qza("woltka_species_taxonomy.qza")$data %>% parse_taxonomy()

taxasums<-summarize_taxa(SVs, taxonomy)$Genus

mouse_barplot<- taxa_barplot(taxasums, metadata, "mouseid") #mouseid is a metadata column 

mouse_barplot


##Supplementary Figure 4B ##

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


# Plot
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

# Pairwise permutation test
read.delim("truseq_deicode_withmd.txt", sep=" ")->truseq

pairwise.perm.t.test(truseq$distance,truseq$group)


## Supplementary Figure 4C ##

#DEICODE implementation of Robust Aitchison Distance calculated in QIIME2

deicode <- read_qza("mice_species_ordination_deicode.qza")


pcoa_mice <- deicode$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) 

pcoa_mice <- as_tibble(pcoa_mice)

pcoa_mice <- ggplot(pcoa_mice, aes(x=PC1, y=PC2, color=`host_subject_id`, group=`host_subject_id`, shape=`library_method`)) +
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
  ggplot(aes(x = PC1, y = PC2, color = host_subject_id, group = host_subject_id, shape = library_method)) +
  geom_mark_ellipse() +
  geom_point(alpha = 0.5) +
  theme_q2r() +
  geom_point(size = 3) +
  scale_color_discrete(name = "host_subject_id") +
  labs(title = "PCoA Plot",
       x = paste0("PC1 (", variance[1], ")"),
       y = paste0("PC2 (", variance[2], ")"))


pcoa_mice


```


