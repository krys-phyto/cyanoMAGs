# Krys Kibler
# Purpose: Pi in genes



# Library 
library(tidyverse)
library(readr)
library(ggpubr)

# Files

names <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
gene.loc <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-locations/gene-calls.txt"
gene.fun <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-coverages/16cyanoMags-refined_FUNCTIONS.txt"

list_of_files <- list.files(path = "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/instrain/profile/",
                            recursive = TRUE,
                            pattern = "*/*gene_info.tsv$",
                            full.names = TRUE)
nuc.diversity <- list_of_files %>%
  set_names() %>% 
  map_df(read_tsv, .id = "file_name")  

# Read in files
names <- read_delim(names,
                    "\t", escape_double = FALSE, trim_ws = TRUE)


nuc.diversity$bins <- gsub('_contig_.+','',nuc.diversity$scaffold)
nuc.diversity$bins <- sub("^", "mag_", nuc.diversity$bins)
nuc.diversity$end <- nuc.diversity %>% with(end + 1)
nuc.diversity$direction <- gsub("-1", "r", nuc.diversity$direction)
nuc.diversity$direction <- gsub("1", "f", nuc.diversity$direction)
nuc.diversity <- left_join(nuc.diversity, gene.loc, by = c("scaffold"="contig", "start", "end"="stop", "direction"))
nuc.diversity <- left_join(nuc.diversity, gene.fun, by = "gene_callers_id")


# Tidying nuc.diversity

nuc.diversity <- left_join(nuc.diversity, names, by = "bins")

nuc.diversity <- nuc.diversity %>% subset(descriptive_nickname != "VAM-591" &
                                            descriptive_nickname != "VAM-558")

nuc.diversity.clean <- nuc.diversity %>% drop_na(nucl_diversity)
nuc.diversity.clean <- nuc.diversity.clean %>% filter(coverage >= 5)

# calc log (pi)
nuc.diversity.clean$log_pi <- log(nuc.diversity.clean$nucl_diversity)


# count
nuc.diversity.clean$count <- 1

nuc.diversity.clean <- nuc.diversity.clean %>% 
  filter(descriptive_nickname == "PLX-501" |
           descriptive_nickname == "DOL-542" |
           descriptive_nickname == "APH-490" |
           descriptive_nickname == "Msyn-517" |
           descriptive_nickname == "Maer-573" |
           descriptive_nickname == "CYM-509") 

nuc.diversity_count <- nuc.diversity.clean %>% 
  group_by(bins, descriptive_nickname) %>% 
  summarise(Freq = sum(count))

# plotting
plot <- nuc.diversity.clean %>% 
  filter(descriptive_nickname == "PLX-501" |
           descriptive_nickname == "DOL-542" |
           descriptive_nickname == "APH-490" |
           descriptive_nickname == "Maer-573") %>% 
  mutate(across(descriptive_nickname, factor, levels=c("PLX-501", "DOL-542", "APH-490", "Maer-573"))) %>%
  ggplot(aes(x = descriptive_nickname, y = nucl_diversity, fill=factor(descriptive_nickname, levels = c("PLX-501", "DOL-542", "APH-490", "Msyn-517", "Maer-573", "CYM-509")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("PLX-501" = "#2980B9",
                               "DOL-542" = "#512E5F",
                               "APH-490" = "#A569BD",
                               "Maer-573" = "#F98E0D",
                               "CYM-509" = "#196F3D")) +
  coord_flip() +
  ylab('Gene Nucleotide Diversity (π)') +
  xlab('Cyano MAGs') +
  labs(fill="Cyano Bins") +
  annotate("text", x= "PLX-501", y=0.06, label= "(10,532)", size = 4) +
  annotate("text", x= "DOL-542", y=0.06, label= "(13,745)", size = 4) +
  annotate("text", x= "APH-490", y=0.06, label= "(20,284)", size = 4) +
  annotate("text", x= "Maer-573", y=0.06, label= "(5646)", size = 4) +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust= 0.5, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none") +
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = "Maer-573", size = 6)  









