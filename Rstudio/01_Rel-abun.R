# Krystyn Kibler 2023-01-11
# Purpose: Plotting relative abundance for paper - Figure 2


### Files ###
file_relabund <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bins_across_samples/abundance.txt"
file_nicknames <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
file_taxonomy <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-taxonomy/SCG-TAXONOMY.txt"

#### Libraries ###
library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(stringr)
library(dplyr)
library(scales)

### Load in data ###
relabund <- read_delim(file_relabund,
                       "\t", escape_double = FALSE, trim_ws = TRUE)
nicknames <- read_delim(file_nicknames,
                        "\t", escape_double = FALSE, trim_ws = TRUE)
scgtax <- read_delim(file_taxonomy,
                     "\t", escape_double = FALSE, trim_ws = TRUE)

# tidy data
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


### Calculations ###
  # relabundance as % 'dominance' of bins
relabund.sum <- relabund %>% group_by(sampledate) %>%
        summarise(relabund_sum = sum(relabund))
relabund <- left_join(relabund, relabund.sum, by = "sampledate")

relabund$perc <- relabund %>% with(relabund / relabund_sum)
relabund$perc <- relabund$perc *100
relabund <- replace(relabund, is.na(relabund), 0)

relabund$sampledate <- as_date(relabund$sampledate)

# year total of relabund

relabund_yeartotal <- relabund %>% group_by(year4, short_name) %>% 
  summarise(sum_perc = sum(perc))
  
year_count <- data.frame(year4 = c(2008,2009,2010,2011,2012),
                         samples = c(16, 15, 21, 13, 27))

relabund_yeartotal <- left_join(relabund_yeartotal, year_count, by = "year4")

relabund_yeartotal$sum_perc_con <- relabund_yeartotal %>% with(sum_perc / samples)

# Tidy relabund.sum
relabund.sum$sampledate <- as_date(relabund.sum$sampledate)
relabund.sum$year4 <- year(relabund.sum$sampledate)
relabund.sum$jday <- yday(relabund.sum$sampledate)
relabund.sum$month <- month(relabund.sum$sampledate)
relabund.sum$day <- as.numeric(day(relabund.sum$sampledate))
relabund.sum$month_day <- as_date(paste(2023,month(relabund.sum$sampledate),day(relabund.sum$sampledate), sep = "-"))

# plots
relabund$date <- gsub("20..-", "", relabund$sampledate)
relabund$date <- gsub("03-", "Mar-", relabund$date)
relabund$date <- gsub("04-", "Apr-", relabund$date)
relabund$date <- gsub("05-", "May-", relabund$date)
relabund$date <- gsub("06-", "Jun-", relabund$date)
relabund$date <- gsub("07-", "Jul-", relabund$date)
relabund$date <- gsub("08-", "Aug-", relabund$date)
relabund$date <- gsub("09-", "Sep-", relabund$date)
relabund$date <- gsub("10-", "Oct-", relabund$date)
relabund$date <- gsub("11-", "Nov-", relabund$date)

dates <- unique(relabund$sampledate)
dates <- gsub("20..-", "", dates)
dates <- gsub("-03-", "Mar-", dates)
dates <- gsub("-04-", "Apr-", dates)
dates <- gsub("-05-", "May-", dates)
dates <- gsub("-06-", "Jun-", dates)
dates <- gsub("-07-", "Jul-", dates)
dates <- gsub("-08-", "Aug-", dates)
dates <- gsub("-09-", "Sep-", dates)
dates <- gsub("-10-", "Oct-", dates)
dates <- gsub("-11-", "Nov-", dates)

custom_labels = c("2008-06-26"="Jun-26", "2008-07-03"="Jul-03", "2008-07-09"="Jul-09", "2008-07-18"="Jul-18", "2008-07-19"="Jul-19",
  "2008-07-21"="Jul-21", "2008-07-23"="Jul-23", "2008-08-05"="Aug-05", "2008-08-13"="Aug-13", "2008-08-20"="Aug-20",
  "2008-08-27"="Aug-27", "2008-09-12"="Sep-12", "2008-09-13"="Sep-13", "2008-09-25"="Sep-25", "2008-10-08"="Oct-08",
  "2008-10-17"="Oct-17", "2009-04-22"="Apr-22", "2009-04-29"="Apr-29", "2009-06-09"="Jun-09", "2009-06-18"="Jun-18",
  "2009-06-26"="Jun-26", "2009-07-07"="Jul-07", "2009-07-30"="Jul-30", "2009-08-10"="Aug-10", "2009-08-26"="Aug-26",
  "2009-09-13"="Sep-13", "2009-09-14"="Sep-14", "2009-09-27"="Sep-27", "2009-10-07"="Oct-07", "2009-10-26"="Oct-26",
  "2009-11-14"="Nov-14", "2010-04-20"="Apr-20", "2010-05-05"="May-05", "2010-05-18"="May-18", "2010-05-20"="May-20",
  "2010-06-02"="Jun-02", "2010-06-13"="Jun-13", "2010-06-15"="Jun-15", "2010-06-21"="Jun-21", "2010-07-06"="Jul-06",
  "2010-07-15"="Jul-15", "2010-07-16"="Jul-16", "2010-07-27"="Jul-27", "2010-08-05"="Aug-05", "2010-08-17"="Aug-17",
  "2010-08-30"="Aug-30", "2010-08-31"="Aug-31", "2010-09-14"="Sep-14", "2010-09-26"="Sep-26", "2010-10-13"="Oct-13",
  "2010-10-29"="Oct-29", "2010-11-19"="Nov-19", "2011-05-03"="May-03", "2011-05-18"="May-18", "2011-06-01"="Jun-01",
  "2011-06-13"="Jun-13", "2011-06-28"="Jun-28", "2011-07-12"="Jul-12", "2011-07-25"="Jul-25", "2011-08-09"="Aug-09",
  "2011-08-22"="Aug-22", "2011-09-04"="Sep-04", "2011-09-21"="Sep-21", "2011-10-03"="Oct-03", "2011-11-01"="Nov-01",
  "2011-11-30"="Nov-30", "2012-03-05"="Apr-05", "2012-04-02"="May-02", "2012-05-05"="May-05", "2012-05-17"="May-17",
  "2012-06-02"="Jun-02", "2012-06-08"="Jun-08", "2012-06-15"="Jun-15", "2012-06-22"="Jun-22", "2012-06-29"="Jun-29",
  "2012-07-06"="Jul-06", "2012-07-13"="Jul-13", "2012-07-17"="Jul-17", "2012-07-20"="Jul-20", "2012-08-03"="Aug-03",
  "2012-08-17"="Aug-17", "2012-08-24"="Aug-24", "2012-08-31"="Aug-31", "2012-09-07"="Sep-07", "2012-09-13"="Sep-13",
  "2012-09-21"="Sep-21", "2012-09-27"="Sep-27", "2012-10-08"="Oct-08", "2012-10-12"="Oct-12", "2012-10-22"="Oct-22",
  "2012-10-26"="Oct-26", "2012-11-05"="Nov-05", "2012-11-09"="Nov-09", "2012-11-16"="Nov-16")

plot1_all <- relabund %>% 
  ggplot(aes(x = as.factor(sampledate), y = perc, fill=factor(descriptive_nickname, levels = c("CYM-577", "CYM-509", "CYM-539", "CYM-551", "CYM-571", "CYM-117", "CYM-118", "VUL-547", "NOD-531", "PLX-501", "DOL-542", "APH-490", "Msyn-517", "Maer-573")))) +
  geom_bar(stat = "identity", color="black", width = 1, size = 0.1) +
  facet_wrap( ~year4 , nrow = 5, scales = "free_x") +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_x_discrete(labels = custom_labels) +
  #scale_x_date(expand = expansion(mult = c(0, 0))) +
  #scale_x_date(date_labels = "%b-%d") +
  labs(y = 'Relative Abundance (Cyano %)', x = 'Date') +
  guides(fill=guide_legend(title="Cyano MAGs")) +
  scale_fill_manual(values = c("CYM-577" = "#145A32",
                               "CYM-509" = "#196F3D",
                               "CYM-539" = "#1E8449",
                               "CYM-551" = "#229954",
                               "CYM-571" = "#27AE60",
                               "CYM-117" = "#52BE80",
                               "CYM-118" = "#7DCEA0",
                               "VUL-547" = "#A9DFBF",
                               "NOD-531" = "#154360",
                               "PLX-501" = "#2980B9",
                               "DOL-542" = "#512E5F",
                               "APH-490" = "#A569BD",
                               "Msyn-517" = "#D35400",
                               "Maer-573" = "#F98E0D")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.text = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_text(face = "bold", size = 18, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = "bold", size = 18),
        legend.position = "left",
        legend.title=element_text(size=18, face = "bold"), 
        legend.text=element_text(size=16)) 
  
#plot_side <- relabund_yeartotal %>% ggplot(aes(x = factor(short_name, levels = c("Microcystaceae", "Nostocaceae", "Phormidesmiaceae", "Cyanobiaceae")), y = sum_perc_con, fill=factor(short_name, levels = c("Cyanobiaceae", "Phormidesmiaceae", "Nostocaceae", "Microcystaceae")))) +
#  geom_bar(stat = "identity", color="black", width = 0.9, size = 0.1) +
#  facet_wrap( ~year4 , nrow = 5) +
#  scale_fill_manual(values = c("Microcystaceae" = "#EC8304",
#                               "Nostocaceae" = "#9C27B0",
#                               "Phormidesmiaceae" = "#1565C0",
#                               "Cyanobiaceae" = "#4CAF50")) +
#  labs(y = 'Sum RA', x = 'Family') +
#  scale_y_continuous(limits = c(0,55), expand = c(0, 0)) +
#  coord_flip() +
#  theme_classic() +
#  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, face = "bold", size = 10),
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.title.y = element_text(size=12, face = "bold"),
#        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
#        axis.title.x = element_text(size=10, face = "bold"),
#        legend.position="none",
#        strip.background = element_blank(),
#        strip.text.x = element_blank())

plot_line <- relabund.sum %>% ggplot(aes(x = month_day, y = relabund_sum)) +
  geom_point() +
  geom_line() +
  labs(x = 'Month', y = 'Relative Abundance (Total %)') +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  scale_x_date(labels = date_format("%b"), breaks = "1 month") +
  facet_wrap( ~year4 , nrow = 5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "bold", size = 14),
        axis.text.y = element_text(angle = 45, hjust=1, face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_text(size=16, face = "bold"),
        axis.title.x = element_text(size=16, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black", face = "bold"))

plot_together <- ggarrange(plot1_all, plot_line,
                           labels = c("A", "B"),
                           ncol = 2, widths = c(5, 1.5),  align = "h")

ggsave("/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/rstudio/figures/fig_1.png", plot = plot_together, width = 6.5, height = 7, units = "in")



# Line plot of cyano relative abundance to cyanos

relabund$sampledate <- as_date(relabund$sampledate)

relabund %>% filter(descriptive_nickname == "APH-490" |
                      descriptive_nickname == "DOL-542" |
                      descriptive_nickname == "Maer-573" |
                      descriptive_nickname == "PLX-501") %>% 
  ggplot(aes(x = sampledate, y = relabund, 
             color = factor(descriptive_nickname, levels = c("APH-490", "DOL-542", "Maer-573", "PLX-501")), 
             group = year(sampledate))) +
  geom_line() +
  geom_point() +
  facet_wrap(~descriptive_nickname, ncol = 1) +
  scale_color_manual(values = c("PLX-501" = "#2980B9",
                                "DOL-542" = "#512E5F",
                                "APH-490" = "#A569BD",
                                "Maer-573" = "#F98E0D")) +
  scale_x_date(date_labels = "%Y",
               breaks = seq.Date(from = as.Date("2008-01-01"), 
                                 to = as.Date("2012-01-01"), 
                                 by = "1 year")) +
  ylab('Relative Abundance (%)') +
  xlab('Year') +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, face = "bold", size = 7),
        axis.text.y = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "none")





coverages <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-coverages/16cyanoMags-refined-COVs.txt"
coverages <- read_delim(coverages)

# Fix headers for columns
coverages <- setNames(coverages, c(metag_header_dates))
names(coverages)[1] <- "contigs"

# Fix contigs names
#582580531_contig_100_split_00001
coverages$bins <- paste("mag_", coverages$contigs, sep = "")
coverages$bins <- gsub("_contig_.+", "", coverages$bins)

# Change original bin names to easy bin names
coverages <- left_join(coverages, nicknames, by = "bins")

# pivot longer
coverages <- coverages %>% pivot_longer(-c(1,96:103), names_to = "sampledate", values_to = "coverage")

# select for those true true cyanos
coverages <- coverages %>% subset(nicknames != "MAG_003" & nicknames != "MAG_012")

coverages$sampledate <- as_date(coverages$sampledate)
coverages$coverage <- as.numeric(coverages$coverage)

coverages.mean <- coverages %>% group_by(bins, descriptive_nickname, sampledate) %>% 
  summarize(coverage = mean(coverage, na.rm = TRUE),
            cov_sd = sd(coverage, na.rm = TRUE))

coverages.mean$count <- 1
count.above5 <- coverages.mean %>% filter(coverage >= 5) %>% 
  summarise(count.above5 = sum(count))
count.above1 <- coverages.mean %>% filter(coverage >= 1) %>% 
  summarise(count.above1 = sum(count))

coverages.mean <- left_join(coverages.mean, count.above5, by = "descriptive_nickname")
coverages.mean <- left_join(coverages.mean, count.above1, by = "descriptive_nickname")

coverages.mean %>% filter(descriptive_nickname == "Maer-573" |
                            descriptive_nickname == "APH-490" |
                            descriptive_nickname == "PLX-501" |
                            descriptive_nickname == "DOL-542") %>% 
  filter(coverage >= 1) %>% 
  ggplot(aes(x = factor(descriptive_nickname, levels = c("Maer-573", "APH-490", "PLX-501", "DOL-542")) , y = coverage, fill = descriptive_nickname)) +
  scale_fill_manual(values = c("PLX-501" = "#2980B9",
                                "DOL-542" = "#512E5F",
                                "APH-490" = "#A569BD",
                                "Maer-573" = "#F98E0D")) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 6, color = "black") +
  geom_text(aes(x = descriptive_nickname, y = 90, label = paste("n =", count.above1, sep = " "), size = 7)) +
  coord_cartesian(ylim = c(0, 90)) +
  labs(x = "Cyano MAGs", y = "Coverage") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.text = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = "bold", size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "none",
        legend.title=element_text(size=16, face = "bold"), 
        legend.text=element_text(size=14)) 

coverages.mean %>% filter(descriptive_nickname == "Maer-573" |
                            descriptive_nickname == "APH-490" |
                            descriptive_nickname == "PLX-501" |
                            descriptive_nickname == "DOL-542") %>% 
  filter(coverage >= 1) %>% group_by(descriptive_nickname) %>% 
  summarise(mean = mean(coverage), 
            sd = sd(coverage))



### Archive ###

# Are there any mags that are consistently below 5% relative abundance
below5 <- relabund %>% subset(perc < 10)
below5$below5 <- ifelse(below5$perc < 5, 1, 0)
below5$below5 <- ifelse(below5$perc < 5, 1, 0)

# These are the bins that across all the metags are 75% consistently below 5% relative abundance
below5_count <- below5 %>% group_by(nicknames) %>%
  summarise(count = sum(below5))

below5_count$class <- ifelse(below5_count$count > 75, "other", "bin")

for (i in 1:nrow(below5_count)) {
  below5_count[i,3] <- ifelse(below5_count[i,3] == "bin", below5_count[i,1], "other")
}

relabund <- left_join(relabund, below5_count, by = "nicknames")

plot2_other <- relabund %>% ggplot(aes(x = as.factor(jday), y = perc, fill = class)) +
  geom_bar(stat = "identity", width = 0.9, size = 2) +
  facet_wrap( ~year4 , nrow = 5) +
  scale_y_continuous(limits = c(0,101), expand = c(0, 0)) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  labs(y = 'Relative Abundance (%)', x = 'Date') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


