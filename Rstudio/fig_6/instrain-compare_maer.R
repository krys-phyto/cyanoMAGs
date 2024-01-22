# Krys Kibler 2023-12-12
# Doin instrain compare things, but on maer-573



library(tidyverse)
library(readr)
library(lubridate)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(plotly)
library(ggpubr)


compare.maer <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/maer-instrain/instrain/instrainComparer/output/instrainComparer_comparisonsTable.tsv"
compare.maer <- read_delim(compare.maer, "\t", escape_double = FALSE, trim_ws = TRUE)

compare.maer.nona <- compare.maer %>% drop_na(conANI)
compare.maer.nona$bins <- gsub("_contig_.+", "", compare.maer.nona$scaffold)



# so things that i need to think about
# percent_genome_compared value, should probably put a threshold on that to keep only values where there is a higher comparison between contigs
compare.avg_2 <- compare.maer.nona %>% filter(percent_genome_compared > .15) %>% 
  group_by(name1, name2) %>% 
  summarize_at(.vars = c("conANI", "popANI", "percent_genome_compared"), .funs = mean)



# pretend that above question does not exist and continue as normal
compare.avg <- compare.maer.nona %>% 
  group_by(name1, name2) %>% 
  summarize_at(.vars = c("conANI", "popANI", "percent_genome_compared"), .funs = mean)

compare.avg$name1_date <- as_date(gsub(".OR.sorted.bam", "", compare.avg$name1))
compare.avg$name2_date <- as_date(gsub(".OR.sorted.bam", "", compare.avg$name2))


plot_compare.all <- compare.avg %>% ggplot(aes(x = as.factor(name1_date), y = as.factor(name2_date), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  xlab("name1 date") +
  ylab("name2 date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(t = 0, r = 15, b = 0, l = 0)))



# RelAbund
relabund <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bins_across_samples/abundance.txt"
relabund <- read_delim(relabund, "\t", escape_double = FALSE, trim_ws = TRUE)

nicknames <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
nicknames <- read_delim(nicknames, "\t", escape_double = FALSE, trim_ws = TRUE)

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
relabund.sum <- relabund %>% 
  group_by(sampledate) %>%
  summarise(relabund_sum = sum(relabund))

relabund <- left_join(relabund, relabund.sum, by = "sampledate")

relabund$perc <- relabund %>% with(relabund / relabund_sum)
relabund$perc <- relabund$perc *100
relabund <- replace(relabund, is.na(relabund), 0)


# Maer Relabund blooms
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

ggplotly(plot_maer.abund)


bloom <- relabund.maer %>% filter(relabund_sum >= 10 & perc >= 70)
# Blooms
# "2008-07-18" "2009-04-22" "2009-09-13" "2009-10-07" "2011-09-04" "2012-06-02"



### Compare with only bloom dates
compare.bloom <- compare.avg %>% filter(name1_date == "2008-07-18" |
                                          name1_date == "2009-04-22" | 
                                          name1_date == "2009-09-13" | 
                                          name1_date == "2010-10-07" | 
                                          name1_date == "2011-09-04" | 
                                          name1_date == "2012-06-02")

compare.bloom <- compare.bloom %>% filter(name2_date == "2008-07-18" |
                                            name2_date == "2009-04-22" | 
                                            name2_date == "2009-09-13" | 
                                            name2_date == "2010-10-07" | 
                                            name2_date == "2011-09-04" | 
                                            name2_date == "2012-06-02")


plot_compare.blooms <- compare.bloom %>% ggplot(aes(x = as.factor(name1_date), y = as.factor(name2_date), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  xlab("Bloom Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        axis.title.y = element_blank())








### bloom date snaps ###

# relative abundance for 2012
relabund.maer$year4 <- year(relabund.maer$sampledate)

plot_maer.abund <- relabund.maer %>% filter(year4 == 2012) %>% 
  ggplot(aes(x = sampledate, y = perc)) +
  geom_line() +
  geom_point() +
  ylab('Perc Maer Dominance (RA)') +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(), 
        legend.position = "none")

ggplotly(plot_maer.abund)

#2012-06-02, 2012-06-29, 2012-07-20, 2012-09-07, 2012-10-12

# new strat for other bloom dates
compare.avg$year4_name1 <- year(compare.avg$name1_date)
compare.avg$year4_name2 <- year(compare.avg$name2_date)
compare.avg$month_name2 <- month(compare.avg$name2_date)
compare.avg$month_name1 <- month(compare.avg$name1_date)


bloom <- relabund.maer %>% filter(sampledate == "2012-06-02" |
                                   sampledate == "2012-06-29" |
                                   sampledate == "2012-07-20" |
                                   sampledate == "2012-09-07" |
                                   sampledate == "2012-10-12")


row.num = nrow(bloom)
datalist = list()

for (i in 1:row.num) {
  
  bloom_date = bloom$sampledate[i]
  bloom_year = year(bloom_date)
  
  compare.bloom <- compare.avg %>% filter(name1_date == bloom_date |
                                            name2_date == bloom_date)
  
  compare.bloom <- compare.bloom %>% filter(year4_name1 == bloom_year)
  compare.bloom <- compare.bloom %>% filter(year4_name2 == bloom_year)
  
  compare.bloom$bloom.compare.date <- bloom_date
  
  row.num.2 = nrow(compare.bloom)
  for (u in 1:row.num.2) { 
    
    date1 = as.character(compare.bloom$name1_date[u])
    date2 = as.character(compare.bloom$name2_date[u])
    compare.bloom$non.bloom.compare.date <- NA
    
    if(compare.bloom$name1_date[u] != bloom_date) compare.bloom$non.bloom.compare.date[u] <- date1 else compare.bloom$non.bloom.compare.date[u] <- date2
    
    compare.bloom$bloom.compare.date <- as_date(compare.bloom$bloom.compare.date)
    compare.bloom$non.bloom.compare.date <- as_date(compare.bloom$non.bloom.compare.date)
  }
  
  compare.bloom[row.num.2 + 1,] <- NA
  compare.bloom$non.bloom.compare.date[row.num.2+1] <- as_date(bloom_date)
  compare.bloom$bloom.compare.date[row.num.2+1] <- as_date(bloom_date)
  
  
  datalist[[i]] <- compare.bloom
}


big.compare <- do.call(rbind, datalist)
rm(datalist)
big.compare <- big.compare[!(is.na(big.compare$name1_date)), ]

row.num = nrow(big.compare)
for (i in 1:row.num) { 
  
  date1 = big.compare$name1_date[i]
  date2 = big.compare$name2_date[i]
  bloom_date = big.compare$bloom.compare.date[i]
  
  
  if (big.compare$name1_date[i] != bloom_date) {
    big.compare$non.bloom.compare.date[i] <- date1
  } else {
    big.compare$non.bloom.compare.date[i] <- date2
  }
  
}

big.compare$bloom.compare.date <- as_date(big.compare$bloom.compare.date)
big.compare$non.bloom.compare.date <- as_date(big.compare$non.bloom.compare.date)

#2012-06-12, 2012-06-29, 2012-07-20, 2012-09-07, 2012-10-12
row.num = nrow(big.compare)
big.compare[row.num + 1,] <- NA
big.compare$non.bloom.compare.date[row.num+1] <- as_date("2012-06-02")
big.compare$bloom.compare.date[row.num+1] <- as_date("2012-06-02")

big.compare[row.num + 2,] <- NA
big.compare$non.bloom.compare.date[row.num+2] <- as_date("2012-06-29")
big.compare$bloom.compare.date[row.num+2] <- as_date("2012-06-29")

big.compare[row.num + 3,] <- NA
big.compare$non.bloom.compare.date[row.num+3] <- as_date("2012-07-20")
big.compare$bloom.compare.date[row.num+3] <- as_date("2012-07-20")

big.compare[row.num + 4,] <- NA
big.compare$non.bloom.compare.date[row.num+4] <- as_date("2012-09-07")
big.compare$bloom.compare.date[row.num+4] <- as_date("2012-09-07")

big.compare[row.num + 5,] <- NA
big.compare$non.bloom.compare.date[row.num+5] <- as_date("2012-10-12")
big.compare$bloom.compare.date[row.num+5] <- as_date("2012-10-12")

big.compare$bloom.compare.date <- as.character(big.compare$bloom.compare.date)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-06-02", 1)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-06-29", 2)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-07-20", 3)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-09-07", 4)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-10-12", 5)

big.compare$bloom.compare.date <- as.numeric(big.compare$bloom.compare.date)


plot_compare.all <- big.compare %>% ggplot(aes(x = non.bloom.compare.date, y = bloom.compare.date, fill = conANI)) +
  geom_hline(yintercept = c(0.5:5.5), linetype = 2, color = 'grey50') +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                      breaks=c(.85,0.925,1.00), 
                      labels=c("0.85","0.925","1.00"),
                      limits=c(.85,1.00))   +
  labs(x = "Month", y = "Peaks") +
  scale_y_continuous(breaks=c(1, 2, 3, 4, 5)) +
  scale_x_date(date_breaks="1 month", date_labels="%b", limits = c(min(big.compare$non.bloom.compare.date), max = max(big.compare$non.bloom.compare.date))) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(face = "bold", size = 14), 
        legend.position = c(0.07, 0.5)) +
  annotate('text', label = 1, x = as.Date('2012-06-02'), y = 5.8, size = 6) +
  annotate('text', label = 2, x = as.Date("2012-06-29"), y = 5.8, size = 6) +
  annotate('text', label = 3, x = as.Date("2012-07-20"), y = 5.8, size = 6) +
  annotate('text', label = 4, x = as.Date("2012-09-07"), y = 5.8, size = 6) +
  annotate('text', label = 5, x = as.Date("2012-10-12"), y = 5.8, size = 6) 

plot_compare.all


plot_together <- ggarrange(plot_maer.abund, plot_compare.all, 
                           ncol = 1, align = "v")



plot_compare.all_scaled <- big.compare %>% ggplot(aes(x = non.bloom.compare.date, y = bloom.compare.date, fill = conANI)) +
  geom_hline(yintercept = c(0.5:5.5), linetype = 2, color = 'grey50') +
  geom_tile() +
  scale_fill_viridis(na.value = "gray30",
                     breaks=c(.92,0.96,1.00), 
                     labels=c("0.92","0.96","1.00"),
                     limits=c(.92,1.00))   +
  labs(x = "Month", y = "Peaks") +
  scale_y_continuous(breaks=c(1, 2, 3, 4, 5)) +
  scale_x_date(date_breaks="1 month", date_labels="%b", limits = c(min(big.compare$non.bloom.compare.date), max = max(big.compare$non.bloom.compare.date))) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text=element_text(size=12, face = "bold"),
        legend.title = element_text(face = "bold", size = 14), 
        legend.position = c(0.05, 0.5)) +
  annotate('text', label = 1, x = as.Date('2012-06-02'), y = 5.8, size = 6) +
  annotate('text', label = 2, x = as.Date("2012-06-29"), y = 5.8, size = 6) +
  annotate('text', label = 3, x = as.Date("2012-07-20"), y = 5.8, size = 6) +
  annotate('text', label = 4, x = as.Date("2012-09-07"), y = 5.8, size = 6) +
  annotate('text', label = 5, x = as.Date("2012-10-12"), y = 5.8, size = 6) 

plot_compare.all_scaled


plot_together <- ggarrange(plot_maer.abund, plot_compare.all_scaled, 
                           ncol = 1, align = "v")
