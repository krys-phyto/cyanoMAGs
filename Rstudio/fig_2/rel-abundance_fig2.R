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

plot1_all <- relabund %>% ggplot(aes(x = as.factor(jday), y = perc, fill=factor(descriptive_nickname, levels = c("CYM-577", "CYM-509", "CYM-539", "CYM-551", "CYM-571", "CYM-117", "CYM-118", "VUL-547", "NOD-531", "PLX-501", "DOL-542", "APH-490", "Msyn-517", "Maer-573")))) +
  geom_bar(stat = "identity", color="black", width = 0.9, size = 0.1) +
  facet_wrap( ~year4 , nrow = 5) +
  scale_y_continuous(limits = c(-1,101), expand = c(0, 0)) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  labs(y = 'Relative Abundance (%)', x = 'Julian Day') +
  guides(fill=guide_legend(title="Cyano Bins")) +
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "left",
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12)) 
  
plot_side <- relabund_yeartotal %>% ggplot(aes(x = factor(short_name, levels = c("Microcystaceae", "Nostocaceae", "Phormidesmiaceae", "Cyanobiaceae")), y = sum_perc_con, fill=factor(short_name, levels = c("Cyanobiaceae", "Phormidesmiaceae", "Nostocaceae", "Microcystaceae")))) +
  geom_bar(stat = "identity", color="black", width = 0.9, size = 0.1) +
  facet_wrap( ~year4 , nrow = 5) +
  scale_fill_manual(values = c("Microcystaceae" = "#EC8304",
                               "Nostocaceae" = "#9C27B0",
                               "Phormidesmiaceae" = "#1565C0",
                               "Cyanobiaceae" = "#4CAF50")) +
  labs(y = 'Sum RA', x = 'Family') +
  scale_y_continuous(limits = c(0,55), expand = c(0, 0)) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 1, face = "bold", size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=12, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        axis.title.x = element_text(size=10, face = "bold"),
        legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

plot_line <- relabund.sum %>% ggplot(aes(x = month_day, y = relabund_sum)) +
  geom_point() +
  geom_line() +
  labs(x = 'Month', y = 'RA (%)') +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  scale_x_date(labels = date_format("%b"), breaks = "1 month") +
  facet_wrap( ~year4 , nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 0.5, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 11),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_text(size=12, face = "bold"),
        axis.title.x = element_text(size=12, face = "bold"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

plot_together <- ggarrange(plot1_all, plot_side, plot_line,
                           labels = c("A", "B", "C"),
                           ncol = 3, widths = c(4, 1, 2),  align = "h")


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


