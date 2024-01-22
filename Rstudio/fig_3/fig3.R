# Krys Kibler
# Purpose: J/Q Nuc diversity through time and reorganize nuc diversity figs



# Library 
library(tidyverse)
library(readr)
library(ggpubr)

# Files

names <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
gene.loc <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-locations/gene-calls.txt"
gene.fun <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-coverages/16cyanoMags-refined_FUNCTIONS.txt"

list_of_files <- list.files(path = "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/all-gene/instrain/profile",
                            recursive = TRUE,
                            pattern = "*/*scaffold_info.tsv$",
                            full.names = TRUE)
nuc.diversity <- list_of_files %>%
  set_names() %>% 
  map_df(read_tsv, .id = "file_name")  


# Read in files
names <- read_delim(names,
                    "\t", escape_double = FALSE, trim_ws = TRUE)

gene.loc <- read_delim(gene.loc,
                       "\t", escape_double = FALSE, trim_ws = TRUE)
gene.loc <- gene.loc %>% select(c(1,2)) 
names(gene.loc)[1] <- "scaffold"

gene.fun <- read_delim(gene.fun,
                       "\t", escape_double = FALSE, trim_ws = TRUE)
gene.fun <- gene.fun %>% filter(source == "COG20_CATEGORY")

nuc.diversity <- left_join(nuc.diversity, gene.loc, by = "scaffold")

# Tidying nuc.diversity
nuc.diversity$bins <- gsub('_contig_.+','',nuc.diversity$contig)
nuc.diversity$bins <- sub("^", "mag_", nuc.diversity$bins)

nuc.diversity <- left_join(nuc.diversity, names, by = "bins")

nuc.diversity <- nuc.diversity %>% subset(descriptive_nickname != "VAM-591" &
                                            descriptive_nickname != "VAM-558")

nuc.diversity <- nuc.diversity %>% drop_na(nucl_diversity)
nuc.diversity <- nuc.diversity %>% filter(coverage_median >= 5)

# calc log (pi)
nuc.diversity$log_pi <- log(nuc.diversity$nucl_diversity)



# count
nuc.diversity$count <- 1

nuc.diversity_count <- nuc.diversity %>% 
  filter(descriptive_nickname == "PLX-501" |
           descriptive_nickname == "DOL-542" |
           descriptive_nickname == "APH-490" |
           descriptive_nickname == "Msyn-517" |
           descriptive_nickname == "Maer-573" |
           descriptive_nickname == "CYM-509") %>% 
  filter(coverage >= 3) %>% 
  filter(log_pi != -Inf)

nuc.diversity_count <- nuc.diversity %>% 
  group_by(bins, descriptive_nickname) %>% 
  summarise(Freq = sum(count))

# plotting
plot <- nuc.diversity %>% 
  filter(descriptive_nickname == "PLX-501" |
           descriptive_nickname == "DOL-542" |
           descriptive_nickname == "APH-490" |
           descriptive_nickname == "Maer-573" |
           descriptive_nickname == "CYM-509") %>% 
  mutate(across(descriptive_nickname, factor, levels=c("Maer-573", "APH-490", "DOL-542", "PLX-501", "CYM-509"))) %>%
  ggplot(aes(x = descriptive_nickname, y = log_pi, fill=factor(descriptive_nickname, levels = c("PLX-501", "DOL-542", "APH-490", "Msyn-517", "Maer-573", "CYM-509")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("PLX-501" = "#2980B9",
                               "DOL-542" = "#512E5F",
                               "APH-490" = "#A569BD",
                               "Maer-573" = "#F98E0D",
                               "CYM-509" = "#196F3D")) +
  coord_flip() +
  ylab(bquote('Log [Pi] (Gene Nuc Diversity)')) +
  xlab('Cyano Bins') +
  labs(fill="Cyano Bins") +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "none") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Maer-573")  


# Annotation with number of gene samples
plot +
  annotate("text", x= "PLX-501", y=-0.5, label= "(2175)") +
  annotate("text", x= "DOL-542", y=-0.5, label= "(2840)") +
  annotate("text", x= "APH-490", y=-0.5, label= "(4425)") +
  annotate("text", x= "Maer-573", y=-0.5, label= "(1253)") 





# nuc.diversity through time # not enough data abandoned

nuc.diversity <- left_join(nuc.diversity, gene.fun, by = c("scaffold" = "gene_callers_id"))


j.q.avg <- nuc.diversity %>% filter(accession == 'Q' |
                                      accession == "J") %>% 
  group_by(accession, file_name, descriptive_nickname) %>% 
  summarize_at(.vars = c("log_pi", "nucl_diversity"), .funs = mean)








