# Krys Kibler
# New figure 5: fuck bgs with antismash

### Libraries ###
library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(stringr)
library(dplyr)
library(scales)

### Files ###

snvs <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-snvs/gene-regions_snv-variability-profile_NUC.txt"
snvs <- read_delim(snvs, "\t", escape_double = FALSE, trim_ws = TRUE)

gene.fun <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-coverages/16cyanoMags-refined_FUNCTIONS.txt"
gene.fun <- read_delim(gene.fun, "\t", escape_double = FALSE, trim_ws = TRUE)

# gene locations
gene.loc <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-locations/gene-calls.txt"
gene.loc <- read_delim(gene.loc, "\t", escape_double = FALSE, trim_ws = TRUE)

# relabund
relabund <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bins_across_samples/abundance.txt"
relabund <- read_delim(relabund, "\t", escape_double = FALSE, trim_ws = TRUE)

# names
nicknames <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
nicknames <- read_delim(nicknames,"\t", escape_double = FALSE, trim_ws = TRUE)


# RelAbund
# current header format = metag_20080719_OR WANT just the date
metag_header_dates <- colnames(relabund)
metag_header_dates <- sapply(strsplit(metag_header_dates, "_"), "[", 2)
metag_header_dates <- as_date(metag_header_dates)

# Fix headers for columns
relabund <- setNames(relabund, c(metag_header_dates))
names(relabund)[1] <- "bins"

# Change original bin names to easy bin names
relabund <- left_join(relabund, nicknames, by = "bins")

# pivot longer
relabund <- relabund %>% pivot_longer(-c(1,96,97,98,99,100,101,102), names_to = "sampledate", values_to = "relabund")

# select for those true true cyanos
relabund <- relabund %>% subset(nicknames != "MAG_003" & nicknames != "MAG_012")


# add year and julianday to relabund
relabund$year4 <- year(relabund$sampledate)
relabund$jday <- yday(relabund$sampledate)

relabund$sampledate <- as_date(relabund$sampledate)


### Calculations ###
# relabundance as % 'dominance' of bins
relabund.sum <- relabund %>% group_by(sampledate) %>%
  summarise(relabund_sum = sum(relabund))
relabund <- left_join(relabund, relabund.sum, by = "sampledate")

relabund$perc <- relabund %>% with(relabund / relabund_sum)
relabund$perc <- relabund$perc *100
relabund <- replace(relabund, is.na(relabund), 0)






### We want Q, J, totalSNVs through time

# calculate number of snvs
snvs$totalSNVs <- snvs %>% with(coverage * departure_from_consensus)

# filter out all the genes with coverage less than 5
snvs <- snvs %>% filter(gene_coverage >= 5)

# genes that had above 5 coverage at any point through time series
#snvs_unique <- snvs 
#snvs_unique$num_samples <- 1
#snvs_unique <- snvs_unique %>% group_by(corresponding_gene_call, contig_name, sample_id, gene_length) %>% 
#  summarise_at(.vars = c("num_samples"), .funs = sum)
#snvs_unique <- left_join(snvs_unique, filter(gene.fun, source == "COG20_CATEGORY"), by = c("corresponding_gene_call" = "gene_callers_id"))
#snvs_unique <- snvs_unique %>% filter(e_value < 0.01)
#write.csv(snvs_unique, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/paper/supp/genes_above5cov.csv", row.names=FALSE)


# check fig 4 is good 
snvs <- left_join(snvs, filter(gene.fun, source == "COG20_CATEGORY"), by = c("corresponding_gene_call" = "gene_callers_id"))

snvs$bins <- gsub("_contig_.+", "", snvs$contig_name)
snvs$bins <- sub("^", "mag_", snvs$bins )

snvs$sample_id <- gsub("metag_", "", snvs$sample_id)
snvs$sample_id <- as_date(snvs$sample_id)

snvs_avg.gene <- snvs %>% group_by(corresponding_gene_call, accession, bins, sample_id) %>% 
  summarise_at(.vars = c("totalSNVs"), .funs = sum)
snvs_avg.gene <- snvs_avg.gene %>% group_by(corresponding_gene_call, accession, bins) %>% 
  summarise_at(.vars = c("totalSNVs"), .funs = mean)

snvs_avg.gene_id <- left_join(snvs_avg.gene, filter(gene.fun, source == "COG20_FUNCTION"), by = c('corresponding_gene_call' = 'gene_callers_id'))
snvs_avg.gene_id <- snvs_avg.gene_id %>% 
  rename(avg_totalSNVs = totalSNVs)
#write.csv(snvs_avg.gene_id, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/paper/paper-1-submission/supp/avgSNVS_gene_ids.csv", row.names=FALSE)


snvs_avg.gene_id <- snvs_avg.gene_id %>% filter(bins == "mag_3300020490_bin6" |
                                                  bins == "mag_3300020573_bin18")

snvs_avg.gene_id <- snvs_avg.gene_id %>% filter(accession.x == 'Q')                                               

snvs_avg.gene_id <- left_join(snvs_avg.gene_id, gene.loc, by = c('corresponding_gene_call' = 'gene_callers_id'))

numbers <- snvs_avg.gene %>% filter(bins == "mag_3300020490_bin6") %>% 
  group_by(accession) %>% 
  summarise_at(.vars = c("totalSNVs"), .funs = mean)
  
plot_subset.gene.regions_aph <- snvs_avg.gene %>% 
  subset(bins == "mag_3300020490_bin6") %>% 
  subset(accession == "P" |
           accession == "O" |
           accession == "M" |
           accession == "L" |
           accession == "H" |
           accession == "E" |
           accession == "C" |
           accession == "V" |
           accession == "X" | #add X
           accession == "Q" |
           accession == "J" ) %>% 
  mutate(accession = factor(accession, levels=c("J", "Q", "V", "X", "C", "E", "H", "L", "M", "O", "P"))) %>%
  ggplot(aes(x = accession, y = totalSNVs, fill=factor(bins, levels = c("mag_3300020490_bin6")))) +
  geom_boxplot() +
 # geom_hline(yintercept = 2.995732) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0,0)) +
  ylab('Average Total SNVs in Genes') +
  xlab('NCBI COG20 Category') +
  labs(fill="Bins") +
  ggtitle("APH-490") +
  scale_x_discrete(labels=str_wrap(c("J" = "(J/SCGs) Translation, ribosomal structure and biogenesis",
                                     "Q" = "(Q/Sec Metabs) Secondary metabolites biosynthesis, transport and catabolism",
                                     "V" = "(V) Defense mechanisms",
                                     "X" = "(X) Mobilome: prophages and transposons",
                                     "C" = "(C) Energy production and conversion",
                                     "E" = "(E) Amino acid transport and metabolism",
                                     "H" = "(H) Coenzyme transport and metabolism",
                                     "L" = "(L) Replication, recombination and repair",
                                     "M" = "(M) Cell wall/membrane/envelope biogenesis",
                                     "O" = "(O) Posttranslational modification, protein turnover, chaperones",
                                     "P" = "(P) Inorganic ion transport and metabolism"), width = 10)) +
  scale_fill_manual(values = c("mag_3300020490_bin6" = "#A569BD")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 8),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.position = "none") 

plot_1 <- plot_subset.gene.regions_aph +
  #stat_compare_means(method = "anova", label.y = 7.5) +
  stat_compare_means(label = "p.signif", method = "t.test", size = 6,
                     ref.group = "J", na.rm = TRUE, label.y = 900)



numbers <- snvs_avg.gene %>% filter(bins == "mag_3300020573_bin18") %>% 
  group_by(accession) %>% 
  summarise_at(.vars = c("totalSNVs"), .funs = mean)

plot_subset.gene.regions_maer <- snvs_avg.gene %>% 
  subset(bins == "mag_3300020573_bin18") %>% 
  subset(accession == "P" |
           accession == "O" |
           accession == "M" |
           accession == "L" |
           accession == "H" |
           accession == "E" |
           accession == "C" |
           accession == "X" |
           accession == "V" |
           accession == "Q" |
           accession == "J" ) %>% 
  mutate(accession = factor(accession, levels=c("J", "Q", "V", "X", "C", "E", "H", "L", "M", "O", "P"))) %>%
  ggplot(aes(x = accession, y = totalSNVs, fill=factor(bins, levels = c("mag_3300020573_bin18")))) +
  geom_boxplot() +
  # geom_hline(yintercept = 2.995732) +
  scale_y_continuous(limits = c(0, 1000), expand = c(0,0)) +
  ylab('Average Total SNVs in Genes') +
  xlab('NCBI COG20 Category') +
  labs(fill="Bins") +
  ggtitle("Maer-573") +
  scale_x_discrete(labels=str_wrap(c("J" = "(J/SCGs)",
                                     "Q" = "(Q/Sec Metabs)",
                                     "V" = "(V)",
                                     "X" = "(X)",
                                     "C" = "(C)",
                                     "E" = "(E)",
                                     "H" = "(H)",
                                     "L" = "(L)",
                                     "M" = "(M)",
                                     "O" = "(O)",
                                     "P" = "(P)"), width = 10)) +
  scale_fill_manual(values = c("mag_3300020573_bin18" = "#F98E0D")) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.position = "none") 

plot_2 <- plot_subset.gene.regions_maer +
  #stat_compare_means(method = "anova", label.y = 7.5) +
  stat_compare_means(label = "p.signif", method = "t.test", size = 6, 
                     ref.group = "J", na.rm = TRUE, label.y = 900)

plot_together <- ggarrange(plot_1, plot_2,
                           nrow = 2, align = "h")

# JQ
snvs_JQ <- snvs %>% 
  filter(accession == "J" |
           accession == "Q")

snvs_JQ_aph <- snvs_JQ %>% filter(bins == "mag_3300020490_bin6")
snvs_JQ_aph$year4 <- year(snvs_JQ_aph$sample_id)
snvs_JQ_aph$totalSNVs_len <- snvs_JQ_aph %>% with(totalSNVs / gene_length)



snvs_JQ_aph_mean <- snvs_JQ_aph %>% group_by(accession, year4, sample_id) %>% 
  summarise_at(.vars = c("totalSNVs_len"), .funs = mean)

  
# total amount of snvs across all genes of the category
snvs_JQ_aph_sum <- snvs_JQ_aph %>% group_by(sample_id, accession, corresponding_gene_call) %>% 
  summarise_at(.vars = c("totalSNVs"), .funs = sum)
snvs_JQ_aph_sum$count <- 1
snvs_JQ_aph_sum <- snvs_JQ_aph_sum %>% group_by(sample_id, accession) %>% 
  summarise_at(.vars = c("totalSNVs", "count"), .funs = sum)

snvs_JQ_aph_sum$totalSNVs_count <- snvs_JQ_aph_sum %>% with(totalSNVs / count)

snvs_JQ_aph_sum$year4 <- year(snvs_JQ_aph_sum$sample_id)

snvs_JQ_aph_sum_2 <- snvs_JQ_aph_sum %>% select(1,2,5,6) %>% pivot_wider(names_from = accession, values_from = totalSNVs_count)
snvs_JQ_aph_sum_2[is.na(snvs_JQ_aph_sum_2)] <- 0

snvs_JQ_aph_sum_2 <- snvs_JQ_aph_sum_2 %>%  pivot_longer(c(3,4), values_to = "totalSNVs_count", names_to = "accession")

snvs_JQ_aph_sum_2 <- full_join(snvs_JQ_aph_sum_2, relabund.aph, by = c("sample_id" = "sampledate"))

snvs_JQ_aph_sum_3 <- snvs_JQ_aph_sum_2 %>% pivot_wider(names_from = accession, values_from = totalSNVs_count) %>% select(1:7)

snvs_JQ_aph_sum_3[is.na(snvs_JQ_aph_sum_3)] <- 0

snvs_JQ_aph_sum_3 <- snvs_JQ_aph_sum_3 %>% pivot_longer(c(6,7), values_to = "totalSNVs_count", names_to = "accession")
snvs_JQ_aph_sum_3$year4 <- year(snvs_JQ_aph_sum_3$sample_id)

# Plotting for current fig 5

plot_snvs <- snvs_JQ_aph_sum_3 %>% filter(year4 == 2012) %>% 
  ggplot(aes(x = sample_id, y = totalSNVs_count, color = accession)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("Q" = "black", "J" = "grey70")) +
  scale_linetype_manual(values = c("Q" = "solid", "J" = "twodash")) +
  ylab("Total SNVs / (number of genes)") +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  scale_y_continuous(limits = c(0,500), expand = c(0, 0)) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        legend.position = c(.1, 0.815),
        legend.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=14, face = "bold"))

  
plot_rel.abund <- relabund %>% filter(descriptive_nickname == "APH-490" &
                                        year4 == 2012) %>% 
  ggplot(aes(x = sampledate, y = relabund)) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  scale_y_continuous(limits = c(0, 12), expand = c(0, 0)) +
  ylab('APH-490 RA (%)') +
  xlab('Month') +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(), 
        legend.text=element_text(size=16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"))


plot_together <- ggarrange(plot_snvs, plot_rel.abund, plot_compare.all, 
                           ncol = 1, align = "v", 
                           font.label = list(size = 16, color = "black"), vjust = -.05)

ggplotly(plot_together)

# dnds

# Load in snv files from instrain (small issue with combining a couple 2008 files - ignore for now, focus on 2012)
list_of_files <- list.files(path = "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/aph-instrain/instrain/profile/",
                            recursive = FALSE,
                            pattern = "2012.+.OR.sam.IS_SNVs.tsv$",
                            full.names = TRUE)
snv.instrain <- list_of_files %>%
  set_names() %>% 
  map_df(read_tsv, .id = "file_name") 



snv.instrain <- snv.instrain %>% filter(mutation_type != "I")
snv.instrain$file_name <- gsub("/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/aph-instrain/instrain/profile//", "", snv.instrain$file_name)
snv.instrain$file_name <- gsub(".OR.sam.IS_SNVs.tsv", "", snv.instrain$file_name)
snv.instrain$file_name <- as_date(snv.instrain$file_name)


snv.instrain <- snv.instrain %>%  full_join(gene.loc, by=c("scaffold" = "contig")) %>%  
  filter(position >= start, position <= stop ) 

snv.instrain <- left_join(snv.instrain, filter(gene.fun, source == "COG20_CATEGORY"), by = "gene_callers_id")

snv.instrain_JQ <- snv.instrain %>% filter(accession == "Q" |
                                          accession == "J")


# calculating dn/ds = (num of n / num of n sites) / (num of s / num of s sites)
snv.instrain_JQ$totalSNVs <- snv.instrain_JQ %>% with(position_coverage * var_freq)
snv.instrain_JQ$count <- 1
snv.instrain_JQ_dnds <- snv.instrain_JQ %>% group_by(file_name, accession, gene_callers_id, mutation_type) %>% 
  summarise_at(.vars = c("totalSNVs", "count"), .funs = sum)

snv.instrain_JQ_dnds$dx <- snv.instrain_JQ_dnds %>% with(totalSNVs / count)

snv.instrain_JQ_dnds <- ungroup(snv.instrain_JQ_dnds)
snv.instrain_JQ_dnds_wide <- snv.instrain_JQ_dnds%>% select(-c(5,6)) %>% 
  pivot_wider(names_from = mutation_type, values_from = dx)



plot_dnds <- snv.instrain_JQ_dnds_wide %>% ggplot(aes(x = S, y = N, color = accession)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  ylab('dN') +
  xlab('dS') +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(size=16, face = "bold"), 
        legend.text=element_text(size=16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"))




### create a data table that has total number of snvs normalized by gene coverage
snvs_avg.gene_norm <- snvs
snvs_avg.gene_norm$snv_positions <- 1

# Sum total snvs across each gene
snvs_avg.gene_norm <- snvs_avg.gene_norm %>% group_by(corresponding_gene_call, accession, bins, sample_id, gene_coverage) %>% 
  summarise_at(.vars = c("totalSNVs", "snv_positions"), .funs = sum)


# average snvs and maintain gene cov
snvs_avg.gene_norm <- snvs_avg.gene_norm %>% group_by(corresponding_gene_call, accession, bins) %>% 
  summarise_at(.vars = c("totalSNVs", "gene_coverage", "snv_positions"), .funs = mean)

# rename columns
snvs_avg.gene_norm <- snvs_avg.gene_norm %>% 
  rename(avg_totalSNVs = totalSNVs)
snvs_avg.gene_norm <- snvs_avg.gene_norm %>% 
  rename(avg_geneCov = gene_coverage)
snvs_avg.gene_norm <- snvs_avg.gene_norm %>% 
  rename(avg_SNVpos = snv_positions)



# normalize totalSNVs by gene coverage and snv positions
snvs_avg.gene_norm$avg_totalSNVs_normGeneCov <- snvs_avg.gene_norm %>% with(avg_totalSNVs / avg_geneCov)
snvs_avg.gene_norm$avg_totalSNVs_normSNVpos <- snvs_avg.gene_norm %>% with(avg_totalSNVs / avg_SNVpos)


# bring back the NCBI COG20 functions
snvs_avg.gene_norm <- left_join(snvs_avg.gene_norm, filter(gene.fun, source == "COG20_FUNCTION"), by = c('corresponding_gene_call' = 'gene_callers_id'))

# save_as_csv
write.csv(snvs_avg.gene_norm, "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/rstudio/fig-4-5_gene-div_sec-metab/avgSNVS_gene_ids_normalized.csv", row.names=FALSE)


