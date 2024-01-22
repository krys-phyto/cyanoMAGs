# Krys Kibler 2023-12-26
# conANI comparisons with instrain for both aph and maer 


### Libraries ###
library(tidyverse)
library(readr)
library(lubridate)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(plotly)
library(ggpubr)



### Files ###

compare.aph <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/aph-instrain/instrain/instrainComparer/output/instrainComparer_comparisonsTable.tsv"
compare.aph <- read_delim(compare.aph,"\t", escape_double = FALSE, trim_ws = TRUE)

compare.maer <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/maer-instrain/instrain/instrainComparer/output/instrainComparer_comparisonsTable.tsv"
compare.maer <- read_delim(compare.maer, "\t", escape_double = FALSE, trim_ws = TRUE)


### Compare Tidying ###
# no nas
compare.aph<- compare.aph %>% drop_na(conANI)
compare.aph$bins <- gsub("_contig_.+", "", compare.aph$scaffold)

compare.maer<- compare.maer %>% drop_na(conANI)
compare.maer$bins <- gsub("_contig_.+", "", compare.maer$scaffold)

# average contig conANIs with eachother
compare.avg.aph <- compare.aph %>%  
  group_by(name1, name2) %>% 
  summarize_at(.vars = c("conANI", "popANI", "percent_genome_compared"), .funs = mean)

compare.avg.maer <- compare.maer %>% 
  group_by(name1, name2) %>% 
  summarize_at(.vars = c("conANI", "popANI", "percent_genome_compared"), .funs = mean)

# name columns to just dates
compare.avg.maer$name1 <- as_date(gsub(".OR.sorted.bam", "", compare.avg.maer$name1))
compare.avg.maer$name2 <- as_date(gsub(".OR.sorted.bam", "", compare.avg.maer$name2))

compare.avg.aph$name1 <- as_date(gsub(".OR.sorted.bam", "", compare.avg.aph$name1))
compare.avg.aph$name2 <- as_date(gsub(".OR.sorted.bam", "", compare.avg.aph$name2))


### Relative Abundance ###
relabund <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bins_across_samples/abundance.txt"
relabund <- read_delim(relabund, "\t", escape_double = FALSE, trim_ws = TRUE)

nicknames <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
nicknames <- read_delim(nicknames,"\t", escape_double = FALSE, trim_ws = TRUE)

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


# relabundance as % 'dominance' of bins
relabund.sum <- relabund %>% 
  group_by(sampledate) %>%
  summarise(relabund_sum = sum(relabund))

relabund <- left_join(relabund, relabund.sum, by = "sampledate")

relabund$perc <- relabund %>% with(relabund / relabund_sum)
relabund$perc <- relabund$perc *100
relabund <- replace(relabund, is.na(relabund), 0)

# APH relative abundance blooms
relabund.aph <- relabund %>% filter(descriptive_nickname == "APH-490")

relabund.aph <- relabund.aph %>% select(c(descriptive_nickname, sampledate, relabund_sum, perc))

relabund.aph_wide <- relabund.aph %>% pivot_longer(c(3,4), names_to = "category", values_to = "relabundance")

relabund.aph_wide$year <- year(relabund.aph_wide$sampledate)

relabund.aph_wide$mday <- format(as.Date(relabund.aph_wide$sampledate), "%m-%d")

plot_aph.abund <- relabund.aph_wide %>% ggplot(aes(x = mday, y = relabundance, color = category)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0)))


# Maer relative abundance blooms
relabund.maer <- relabund %>% filter(descriptive_nickname == "Maer-573")

relabund.maer <- relabund.maer %>% select(c(descriptive_nickname, sampledate, relabund_sum, perc))

relabund.maer_wide <- relabund.maer %>% pivot_longer(c(3,4), names_to = "category", values_to = "relabundance")

relabund.maer_wide$year <- year(relabund.maer_wide$sampledate)

relabund.maer_wide$mday <- format(as.Date(relabund.maer_wide$sampledate), "%m-%d")

plot_maer.abund <- relabund.maer_wide %>% ggplot(aes(x = mday, y = relabundance, color = category)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0)))


# interactive relative abundance plots
ggplotly(plot_aph.abund)
ggplotly(plot_maer.abund)




### Identify blooms ###
blooms.aph <- relabund.aph %>% filter(relabund_sum >= 10 & perc >= 70)
blooms.maer <- relabund.maer %>% filter(relabund_sum >= 10 & perc >= 70)


compare.avg.aph$year4_name1 <- year(compare.avg.aph$name1)
compare.avg.aph$year4_name2 <- year(compare.avg.aph$name2)
compare.avg.aph$month_name2 <- month(compare.avg.aph$name2)
compare.avg.aph$month_name1 <- month(compare.avg.aph$name1)

compare.avg.maer$year4_name1 <- year(compare.avg.maer$name1)
compare.avg.maer$year4_name2 <- year(compare.avg.maer$name2)
compare.avg.maer$month_name2 <- month(compare.avg.maer$name2)
compare.avg.maer$month_name1 <- month(compare.avg.maer$name1)


# bloom conANI
blooms.aph <- left_join(blooms.aph, compare.avg.aph, by = c("sampledate" = "name1"))

blooms.aph <- blooms.aph %>% filter(name2 == "2008-07-23" |
                                      name2 == "2008-08-05" |
                                      name2 == "2008-09-13" |
                                      name2 == "2008-10-17" |
                                      name2 == "2009-06-09" |
                                      name2 == "2009-10-26" |
                                      name2 == "2010-05-05" |
                                      name2 == "2010-06-13" |
                                      name2 == "2011-05-03" |
                                      name2 == "2011-06-13" |
                                      name2 == "2012-07-06" |
                                      name2 == "2012-09-13" )


plot_blooms.aph <- blooms.aph %>% ggplot(aes(x = as.factor(sampledate), y = as.factor(name2), fill = conANI)) +
  geom_tile() +
  ggtitle("APH-490 Bloom conANI") +
  scale_fill_viridis(na.value = "gray30",
                     breaks=c(.97,0.985,1.00), 
                     labels=c("0.97","0.985","1.00"),
                     limits=c(.97,1.00)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        axis.ticks = element_blank())

#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
#axis.text.y = element_text(face = "bold", size = 12)

#axis.text.x = element_blank(),
#axis.text.y = element_blank(),
#axis.ticks = element_blank()

blooms.maer <- left_join(blooms.maer, compare.avg.maer, by = c("sampledate" = "name1"))

blooms.maer <- blooms.maer %>% filter(name2 == "2008-07-18" |
                                      name2 == "2009-04-22" |
                                      name2 == "2009-09-13" |
                                      name2 == "2009-10-07" |
                                      name2 == "2011-09-04" |
                                      name2 == "2012-06-02" )

# date switching

blooms.maer$sampledate[12] <- as_date("2009-04-22") 
blooms.maer$sampledate[13] <- as_date("2009-09-13")
blooms.maer$sampledate[14] <- as_date("2009-10-07")
blooms.maer$sampledate[15] <- as_date("2011-09-04")

blooms.maer$name2[12] <- as_date("2012-06-02")
blooms.maer$name2[13] <- as_date("2012-06-02")
blooms.maer$name2[14] <- as_date("2012-06-02")
blooms.maer$name2[15] <- as_date("2012-06-02")
                                      
plot_blooms.maer <- blooms.maer %>% ggplot(aes(x = as.factor(sampledate), y = as.factor(name2), fill = conANI)) +
  geom_tile() +
  ggtitle("Maer-573 Bloom conANI") +
  scale_fill_viridis(na.value = "gray30",
                     breaks=c(.97,0.985,1.00), 
                     labels=c("0.97","0.985","1.00"),
                     limits=c(.97,1.00)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        axis.ticks = element_blank())


# plot together

ggarrange(plot_blooms.aph, plot_blooms.maer, ncol = 2, common.legend = TRUE, legend = 'right')


### Identify shoulders 
non.blooms.aph <- relabund.aph %>% filter(perc <= 70)
non.blooms.maer <- relabund.maer %>% filter(perc <= 70)


# nonbloom conANI
non.blooms.aph <- left_join(non.blooms.aph, compare.avg.aph, by = c("sampledate" = "name1"))

non.blooms.aph <- non.blooms.aph %>% filter(name2 != "2008-07-23" |
                                      name2 != "2008-08-05" |
                                      name2 != "2008-09-13" |
                                      name2 != "2008-10-17" |
                                      name2 != "2009-06-09" |
                                      name2 != "2009-10-26" |
                                      name2 != "2010-05-05" |
                                      name2 != "2010-06-13" |
                                      name2 != "2011-05-03" |
                                      name2 != "2011-06-13" |
                                      name2 != "2012-07-06" |
                                      name2 != "2012-09-13" )


plot_nonblooms.aph <- non.blooms.aph %>% ggplot(aes(x = as.factor(sampledate), y = as.factor(name2), fill = conANI)) +
  geom_tile() +
  ggtitle("APH-490 NonBloom conANI") +
  scale_fill_viridis(na.value = "gray30",
                   breaks=c(.60,0.80,1.00), 
                   labels=c("0.60","0.80","1.00"),
                   limits=c(0.60,1.00)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        axis.ticks = element_blank())

#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
#axis.text.y = element_text(face = "bold", size = 12)

#axis.text.x = element_blank(),
#axis.text.y = element_blank(),
#axis.ticks = element_blank()

non.blooms.maer <- left_join(non.blooms.maer, compare.avg.maer, by = c("sampledate" = "name1"))

non.blooms.maer <- non.blooms.maer %>% filter(name2 != "2008-07-18" |
                                        name2 != "2009-04-22" |
                                        name2 != "2009-09-13" |
                                        name2 != "2009-10-07" |
                                        name2 != "2011-09-04" |
                                        name2 != "2012-06-02" )

plot_nonblooms.maer <- non.blooms.maer %>% ggplot(aes(x = as.factor(sampledate), y = as.factor(name2), fill = conANI)) +
  geom_tile() +
  ggtitle("Maer-573 NonBloom conANI") +
  scale_fill_viridis(na.value = "gray30",
                     breaks=c(.60,0.80,1.00), 
                     labels=c("0.60","0.80","1.00"),
                     limits=c(.60,1.00)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        axis.ticks = element_blank())

