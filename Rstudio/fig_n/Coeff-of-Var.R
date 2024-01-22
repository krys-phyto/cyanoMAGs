# Krys Kibler
# Coefficient of Variances for gene coverage for Maer and Aph within sample dates
# CV =  σ / μ
# Coefficient of Variance = std deviation / mean 

# gene functions
gene.fun <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-coverages/16cyanoMags-refined_FUNCTIONS.txt"
gene.fun <- read_delim(gene.fun, "\t", escape_double = FALSE, trim_ws = TRUE)

# gene locations
gene.loc <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-locations/gene-calls.txt"
gene.loc <- read_delim(gene.loc, "\t", escape_double = FALSE, trim_ws = TRUE)

# names
nicknames <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
nicknames <- read_delim(nicknames,"\t", escape_double = FALSE, trim_ws = TRUE)

# gene coverages
gene.cov <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-coverages/16cyanoMags-refined-GENE-COVERAGES.txt"
gene.cov <- read_delim(gene.cov,"\t", escape_double = FALSE, trim_ws = TRUE)


### Libraries ###
library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(stringr)
library(dplyr)
library(scales)


### Tidying ###

# current header format = metag_20080719_OR WANT just the date
genecov_header_dates <- colnames(gene.cov)
genecov_header_dates <- sapply(strsplit(genecov_header_dates, "_"), "[", 2)
genecov_header_dates <- as_date(genecov_header_dates)

# Fix headers for columns
gene.cov <- setNames(gene.cov, c(genecov_header_dates))
names(gene.cov)[1] <- "gene_callers_id"

# pivot longer
gene.cov <- gene.cov %>% pivot_longer(-c(1), names_to = "sample_id", values_to = "gene.cov")

gene.cov$sample_id <- as_date(gene.cov$sample_id)

# pick out genes that are only aph and maer
gene.cov <- left_join(gene.cov, select(gene.loc, c(1,2)), by = "gene_callers_id")

# Bin names and better names
gene.cov$bins <- gsub("_contig_.+", "", gene.cov$contig)
gene.cov$bins <- sub("^", "mag_", gene.cov$bins )

gene.cov <- left_join(gene.cov, nicknames, by = "bins")

gene.cov_maer.aph <- gene.cov %>% filter(descriptive_nickname == "APH-490" |
                                           descriptive_nickname == "Maer-573")

CV_bygene_maer.aph <- gene.cov_maer.aph %>% group_by(descriptive_nickname, gene_callers_id) %>% 
  summarise(sd_gene.cov = sd(gene.cov, na.rm = TRUE),
            mean_gene.cov = mean(gene.cov, na.rm = TRUE))
CV_bygene_maer.aph$CV <- CV_bygene_maer.aph %>% with(sd_gene.cov / mean_gene.cov)

CV_bygene_maer.aph <- left_join(CV_bygene_maer.aph, filter(gene.fun, source == "COG20_FUNCTION"), by = 'gene_callers_id')
write.csv(CV_maer.aph, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/paper/paper-1-submission/supp/CV_mear-aph_bygene.csv", row.names=FALSE)


### Calculations ###
CV_maer.aph <- gene.cov_maer.aph %>% group_by(descriptive_nickname, sample_id) %>% 
  summarise(sd_gene.cov = sd(gene.cov, na.rm = TRUE),
            mean_gene.cov = mean(gene.cov, na.rm = TRUE))

CV_maer.aph$CV <- CV_maer.aph %>% with(sd_gene.cov / mean_gene.cov)

write.csv(CV_maer.aph, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/paper/paper-1-submission/supp/CV_mear-aph.csv", row.names=FALSE)


plot <- CV_maer.aph %>% filter(mean_gene.cov > 5) %>% 
  ggplot(aes(x = as.factor(descriptive_nickname), y = CV, fill = descriptive_nickname)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Maer-573" = "#F98E0D", 
                               "APH-490" = "#A569BD")) +
  ylab('CV (Gene Cov within Metagenomes)') +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(), 
        legend.position = "none")

plot +
  #stat_compare_means(method = "anova", label.y = 7.5) +
  stat_compare_means(method = "t.test", size = 6,
                     na.rm = TRUE, label.y = 1.25)



