# Krys Kibler 2023-07-24
# Doin instrain compare things

library(tidyverse)
library(readr)
library(lubridate)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(plotly)
library(ggpubr)



# Aph-490 and Maer-573 conani/popani before and after bloom 

compare <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/aph-instrain/instrain/instrainComparer/output/instrainComparer_comparisonsTable.tsv"

compare <- read_delim(compare,
                      "\t", escape_double = FALSE, trim_ws = TRUE)

compare.nona <- compare %>% drop_na(conANI)
compare.nona$bins <- gsub("_contig_.+", "", compare.nona$scaffold)


compare.avg <- compare.nona %>% 
  group_by(name1, name2) %>% 
  summarize_at(.vars = c("conANI", "popANI"), .funs = mean)

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

relabund <- read_delim(relabund,
                       "\t", escape_double = FALSE, trim_ws = TRUE)

nicknames <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"

nicknames <- read_delim(nicknames,
                        "\t", escape_double = FALSE, trim_ws = TRUE)

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


# APH Relabund blooms
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

ggplotly(plot_aph.abund)


bloom <- relabund.aph %>% filter(relabund_sum >= 10 & perc >= 70)
# Blooms
#"2008-07-23", "2008-08-05", "2008-09-13", "2008-10-17",
#"2009-06-09", "2009-10-26", "2010-05-05", "2010-06-13"
#"2011-05-03", "2011-06-13", "2012-07-06", "2012-09-13"



### Compare with only bloom dates
compare.bloom <- compare.avg %>% filter(name1_date == "2008-07-23" |
                                    name1_date == "2008-08-05" |
                                    name1_date == "2008-09-13" | 
                                    name1_date == "2008-10-17" |
                                    name1_date == "2009-06-09" | 
                                    name1_date == "2009-10-26" | 
                                    name1_date == "2010-05-05" | 
                                    name1_date == "2010-06-13" |
                                    name1_date == "2011-05-03" | 
                                    name1_date == "2011-06-13" | 
                                    name1_date == "2012-07-06" | 
                                    name1_date == "2012-09-13")

compare.bloom <- compare.bloom %>% filter(name2_date == "2008-07-23" |
                                      name2_date == "2008-08-05" |
                                      name2_date == "2008-09-13" | 
                                      name2_date == "2008-10-17" |
                                      name2_date == "2009-06-09" | 
                                      name2_date == "2009-10-26" | 
                                      name2_date == "2010-05-05" | 
                                      name2_date == "2010-06-13" |
                                      name2_date == "2011-05-03" | 
                                      name2_date == "2011-06-13" | 
                                      name2_date == "2012-07-06" | 
                                      name2_date == "2012-09-13")


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
# new strat for other bloom dates
compare.avg$year4_name1 <- year(compare.avg$name1_date)
compare.avg$year4_name2 <- year(compare.avg$name2_date)
compare.avg$month_name2 <- month(compare.avg$name2_date)
compare.avg$month_name1 <- month(compare.avg$name1_date)


bloom <- relabund.aph %>% filter(sampledate == "2012-06-08" |
                                   sampledate == "2012-07-06" |
                                   sampledate == "2012-07-20" |
                                   sampledate == "2012-09-13" |
                                   sampledate == "2012-10-26")


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
  
 # plot_compare.ind.bloom <- compare.bloom %>% ggplot(aes(x = as.factor(non.bloom.compare.date), y = as.factor(bloom.compare.date), fill = conANI)) +
 #   geom_tile() +
  #  scale_fill_viridis(discrete=FALSE) +
   # labs(x = paste("Bloom", bloom_date, sep = " ")) +
    #theme_classic() +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 0.5, face = "bold", size = 12),
     #     axis.text.y = element_text(face = "bold", size = 12),
      #    panel.border = element_rect(color = "black", fill = NA, size = 1),
       #   axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 12),
        #  axis.title.y = element_blank())

  #ggsave(plot_compare.ind.bloom, file=paste0("/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/rstudio/fig-7_compare/plot_bloom_", bloom_date,".png"), width = 14, height = 10, units = "cm")
  
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

row.num = nrow(big.compare)
big.compare[row.num + 1,] <- NA
big.compare$non.bloom.compare.date[row.num+1] <- as_date("2012-06-08")
big.compare$bloom.compare.date[row.num+1] <- as_date("2012-06-08")

big.compare[row.num + 2,] <- NA
big.compare$non.bloom.compare.date[row.num+2] <- as_date("2012-07-06")
big.compare$bloom.compare.date[row.num+2] <- as_date("2012-07-06")

big.compare[row.num + 3,] <- NA
big.compare$non.bloom.compare.date[row.num+3] <- as_date("2012-07-20")
big.compare$bloom.compare.date[row.num+3] <- as_date("2012-07-20")

big.compare[row.num + 4,] <- NA
big.compare$non.bloom.compare.date[row.num+4] <- as_date("2012-09-13")
big.compare$bloom.compare.date[row.num+4] <- as_date("2012-09-13")

big.compare[row.num + 5,] <- NA
big.compare$non.bloom.compare.date[row.num+5] <- as_date("2012-10-26")
big.compare$bloom.compare.date[row.num+5] <- as_date("2012-10-26")

big.compare$bloom.compare.date <- as.character(big.compare$bloom.compare.date)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-06-08", 1)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-07-06", 2)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-07-20", 3)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-09-13", 4)
big.compare$bloom.compare.date <- replace(big.compare$bloom.compare.date, big.compare$bloom.compare.date == "2012-10-26", 5)

big.compare$bloom.compare.date <- as.numeric(big.compare$bloom.compare.date)

plot_compare.all <- big.compare %>% ggplot(aes(x = non.bloom.compare.date, y = bloom.compare.date, fill = conANI)) +
  geom_hline(yintercept = c(0.5:5.5), linetype = 2, color = 'grey50') +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
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
   annotate('text', label = 1, x = as.Date('2012-06-08'), y = 5.8, size = 6) +
   annotate('text', label = 2, x = as.Date("2012-07-06"), y = 5.8, size = 6) +
   annotate('text', label = 3, x = as.Date("2012-07-20"), y = 5.8, size = 6) +
   annotate('text', label = 4, x = as.Date("2012-09-13"), y = 5.8, size = 6) +
   annotate('text', label = 5, x = as.Date("2012-10-26"), y = 5.8, size = 6) 

plot_compare.all

ggplotly(plot_compare.all)












big.compare$year4 <- year(big.compare$bloom.compare.date)

plot_compare.bloom <- big.compare %>% 
  ggplot(aes(x = as.factor(non.bloom.compare.date), y = as.factor(bloom.compare.date), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "transparent",
                       breaks=c(.92,0.96,1.00), 
                       labels=c("0.92","0.96","1.00"),
                       limits=c(.92,1.00))  +
  labs(y = "Bloom Dates", x = "Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust= 1, face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"))

compare.legend <- get_legend(plot_compare.bloom)


plot_compare.bloom_2008 <- big.compare %>% filter(year4 == 2008) %>% 
  ggplot(aes(x = as.factor(non.bloom.compare.date), y = factor(bloom.compare.date, levels = c("2008-10-17", "2008-09-13", "2008-08-05", "2008-07-23")), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.92,0.96,1.00), 
                     labels=c("0.92","0.96","1.00"),
                     limits=c(.92,1.00))  +
  labs(y = "Bloom Dates", x = "Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust= 1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.key.size = unit(1, 'cm'))

plot_compare.bloom_2009 <- big.compare %>% filter(year4 == 2009) %>% 
  ggplot(aes(x = as.factor(non.bloom.compare.date), y = factor(bloom.compare.date, levels = c("2009-10-26", "2009-06-09")), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.92,0.96,1.00), 
                     labels=c("0.92","0.96","1.00"),
                     limits=c(.92,1.00))  +
  labs(y = "Bloom Dates", x = "Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust= 1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"))

plot_compare.bloom_2010 <- big.compare %>% filter(year4 == 2010) %>% 
  ggplot(aes(x = as.factor(non.bloom.compare.date), y = factor(bloom.compare.date, levels = c("2010-06-13", "2010-05-05")), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.92,0.96,1.00), 
                     labels=c("0.92","0.96","1.00"),
                     limits=c(.92,1.00))  +
  labs(y = "Bloom Dates", x = "Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust= 1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size=16, face = "bold"),
        legend.text=element_text(size=16, face = "bold"))

plot_compare.bloom_2011 <- big.compare %>% filter(year4 == 2011) %>% 
  ggplot(aes(x = as.factor(non.bloom.compare.date), y = factor(bloom.compare.date, levels = c("2011-06-13", "2011-05-03")), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.92,0.96,1.00), 
                     labels=c("0.92","0.96","1.00"),
                     limits=c(.92,1.00))  +
  labs(y = "Bloom Dates", x = "Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust= 1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"))

plot_compare.bloom_2012 <- big.compare %>% filter(year4 == 2012) %>% 
  ggplot(aes(x = as.factor(non.bloom.compare.date), y = factor(bloom.compare.date, levels = c("2012-09-13", "2012-07-06")), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.92,0.96,1.00), 
                     labels=c("0.92","0.96","1.00"),
                     limits=c(.92,1.00))  +
  labs(y = "Bloom Dates", x = "All Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust= 1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 16),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"))

plot.all.years <- ggarrange(plot_compare.bloom_2008, plot_compare.bloom_2009, plot_compare.bloom_2010, 
                            plot_compare.bloom_2011, plot_compare.bloom_2012, nrow = 5, common.legend=TRUE, 
                            legend = "right") 
  



# quick shot for new figure 6 plot
bloom <- c("2012-06-02", "2012-06-08", "2012-06-15", "2012-06-22",
           "2012-06-29", "2012-07-06", "2012-07-13", "2012-07-20",
           "2012-08-03", "2012-08-17", "2012-08-31", "2012-09-13", 
           "2012-09-21", "2012-10-12", "2012-10-26", "2012-11-05")
bloom <- as.data.frame(bloom)
bloom$bloom <- as_date(bloom$bloom)

relabund.aph$year4 <- year(relabund.aph$sampledate)
bloom <- relabund.aph %>% filter(year4 == 2012) %>% filter(sampledate == "2012-06-08" |
                                                             sampledate == "2012-07-06" |
                                                             sampledate == "2012-07-20" |
                                                             sampledate == "2012-10-12")


compare.set <- compare.avg %>% filter(name1_date == "2012-06-08" |
                                        name1_date == "2012-07-06" |
                                        name1_date == "2012-07-20" |
                                        name1_date == "2012-09-13" |
                                        name1_date == "2012-10-26" |
                                        name2_date == "2012-06-08" |
                                        name2_date == "2012-07-06" |
                                        name2_date == "2012-07-20" |
                                        name2_date == "2012-09-13" |
                                        name2_date == "2012-10-26" ) %>% 
                                filter(year4_name1 == 2012)

compare.set.0608 <- compare.set %>% filter(name1_date == "2012-06-08" |
                                             name2_date == "2012-06-08") %>% 
                                    filter(month_name1 == 06 &
                                             month_name2 == 06)
compare.set.0608$bloom.date <- as_date("2012-06-08")
compare.set.0608$non.bloom.date <- NA
compare.set.0608$non.bloom.date[1] <- "2012-06-02"
compare.set.0608$non.bloom.date[2] <- "2012-06-15"
compare.set.0608$non.bloom.date[3] <- "2012-06-29"

nrow = nrow(compare.set.0608)

compare.set.0608[nrow + 1,] <- NA
compare.set.0608$non.bloom.date[nrow + 1] <- "2012-06-08"
compare.set.0608$bloom.date[nrow +1] <- "2012-06-08"
compare.set.0608$non.bloom.date <- as_date(compare.set.0608$non.bloom.date)


plot0608 <- compare.set.0608 %>%  
  ggplot(aes(x = as.factor(non.bloom.date), y = bloom.date, fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.96,0.98,1.00), 
                     labels=c("0.96","0.98","1.00"),
                     limits=c(.96,1.00))  +
  scale_x_discrete(labels=c('Jun-6', 'Jun-8', 'Jun-15', 'Jun-29')) +
  labs(y = "Bloom Dates", x = "All Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(2, 0.5, 2, 0.5), 
                           "inches"),
        axis.ticks.y = element_blank())

compare.set.0706 <- compare.set %>% filter(name1_date == "2012-07-06" |
                                             name2_date == "2012-07-06") %>% 
                                    filter(name1_date == "2012-06-22" |
                                             name1_date == "2012-06-29" |
                                           name2_date == "2012-07-13")
compare.set.0706$bloom.date <- as_date("2012-07-06")
compare.set.0706$non.bloom.date <- NA
compare.set.0706$non.bloom.date[1] <- "2012-06-22"
compare.set.0706$non.bloom.date[2] <- "2012-06-29"
compare.set.0706$non.bloom.date[3] <- "2012-07-13"

nrow = nrow(compare.set.0706)

compare.set.0706[nrow + 1,] <- NA
compare.set.0706$non.bloom.date[nrow + 1] <- "2012-07-06"
compare.set.0706$bloom.date[nrow +1] <- "2012-07-06"
compare.set.0706$non.bloom.date <- as_date(compare.set.0706$non.bloom.date)


plot0706 <- compare.set.0706 %>%  
  ggplot(aes(x = as.factor(non.bloom.date), y = bloom.date, fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.96,0.98,1.00), 
                     labels=c("0.96","0.98","1.00"),
                     limits=c(.96,1.00))  +
  scale_x_discrete(labels=c('Jun-22', 'Jun-29', 'Jul-6', 'Jul-13')) +
  labs(y = "Bloom Dates", x = "All Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(2, 0.5, 2, 0.5), 
                           "inches"),
        axis.ticks.y = element_blank())


compare.set.0720 <- compare.set %>% filter(name1_date == "2012-07-20" |
                                             name2_date == "2012-07-20") %>% 
                                    filter(name1_date == "2012-07-13" |
                                             name2_date == "2012-08-03")
compare.set.0720$bloom.date <- as_date("2012-07-20")
compare.set.0720$non.bloom.date <- NA
compare.set.0720$non.bloom.date[1] <- "2012-07-13"
compare.set.0720$non.bloom.date[2] <- "2012-08-03"

nrow = nrow(compare.set.0720)

compare.set.0720[nrow + 1,] <- NA
compare.set.0720$non.bloom.date[nrow + 1] <- "2012-07-20"
compare.set.0720$bloom.date[nrow +1] <- "2012-07-20"
compare.set.0720$non.bloom.date <- as_date(compare.set.0720$non.bloom.date)


plot0720 <- compare.set.0720 %>%  
  ggplot(aes(x = as.factor(non.bloom.date), y = bloom.date, fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.96,0.98,1.00), 
                     labels=c("0.96","0.98","1.00"),
                     limits=c(.96,1.00))  +
  scale_x_discrete(labels=c('Jul-13', 'Jul-20', 'Aug-3')) +
  labs(y = "Bloom Dates", x = "All Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(2, 0.5, 2, 0.5), 
                           "inches"),
        axis.ticks.y = element_blank())


compare.set.0913 <- compare.set %>% filter(name1_date == "2012-09-13" |
                                             name2_date == "2012-09-13") %>%
                                    filter(name1_date == "2012-08-03" |
                                             name1_date == "2012-08-17" |
                                             name1_date == "2012-08-31" |
                                             name2_date == "2012-09-21")
compare.set.0913$bloom.date <- as_date("2012-09-13")
compare.set.0913$non.bloom.date <- NA
compare.set.0913$non.bloom.date[1] <- "2012-08-13"
compare.set.0913$non.bloom.date[2] <- "2012-08-17"
compare.set.0913$non.bloom.date[3] <- "2012-08-31"
compare.set.0913$non.bloom.date[4] <- "2012-09-21"

nrow = nrow(compare.set.0913)

compare.set.0913[nrow + 1,] <- NA
compare.set.0913$non.bloom.date[nrow + 1] <- "2012-09-13"
compare.set.0913$bloom.date[nrow +1] <- "2012-09-13"
compare.set.0913$non.bloom.date <- as_date(compare.set.0913$non.bloom.date)


plot0913 <- compare.set.0913 %>%  
  ggplot(aes(x = as.factor(non.bloom.date), y = bloom.date, fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.96,0.98,1.00), 
                     labels=c("0.96","0.98","1.00"),
                     limits=c(.96,1.00))  +
  scale_x_discrete(labels=c('Aug-13', 'Aug-17', 'Aug-31', 'Sep-13', 'Sep-21')) +
  labs(y = "Bloom Dates", x = "All Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(2, 0.5, 2, 0.5), 
                           "inches"),
        axis.ticks.y = element_blank())

compare.set.1026 <- compare.set %>% filter(name1_date == "2012-10-26" |
                                             name2_date == "2012-10-26") %>%
                                    filter(name1_date == "2012-10-12" |
                                             name2_date == "2012-11-05")
compare.set.1026$bloom.date <- as_date("2012-10-26")
compare.set.1026$non.bloom.date <- NA
compare.set.1026$non.bloom.date[1] <- "2012-10-12"
compare.set.1026$non.bloom.date[2] <- "2012-11-05"

nrow = nrow(compare.set.1026)

compare.set.1026[nrow + 1,] <- NA
compare.set.1026$non.bloom.date[nrow + 1] <- "2012-10-26"
compare.set.1026$bloom.date[nrow +1] <- "2012-10-26"
compare.set.1026$non.bloom.date <- as_date(compare.set.1026$non.bloom.date)


plot1026 <- compare.set.1026 %>%  
  ggplot(aes(x = as.factor(non.bloom.date), y = bloom.date, fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.96,0.98,1.00), 
                     labels=c("0.96","0.98","1.00"),
                     limits=c(.96,1.00))  +
  scale_x_discrete(labels=c('Oct-12', 'Oct-26', 'Nov-5')) +
  labs(y = "Bloom Dates", x = "All Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(2, 0.5, 2, 0.5), 
                           "inches"),
        axis.ticks.y = element_blank())

compare.set.peaks <- compare.set %>% filter(name1_date != "2012-03-05" &
                                              name1_date != "2012-04-02" &
                                              name1_date != "2012-06-02" &
                                              name1_date != "2012-06-15" &
                                              name1_date != "2012-06-22" &
                                              name1_date != "2012-06-29" &
                                              name1_date != "2012-07-13" &
                                              month_name1 != 8 &
                                              name1_date != "2012-09-21" &
                                              name1_date != "2012-10-08" &
                                              name1_date != "2012-10-12" &
                                              name2_date != "2012-03-05" &
                                              name2_date != "2012-04-02" &
                                              name2_date != "2012-06-02" &
                                              name2_date != "2012-06-15" &
                                              name2_date != "2012-06-22" &
                                              name2_date != "2012-06-29" &
                                              name2_date != "2012-07-13" &
                                              month_name2 != 8 &
                                              name2_date != "2012-09-21" &
                                              name2_date != "2012-10-08" &
                                              name2_date != "2012-10-12" &
                                              month_name2 != 11)

compare.set.all <- rbind(compare.set.0608, compare.set.0706, compare.set.0720, 
                         compare.set.0913, compare.set.1026)

plot1 <- compare.set.all %>%  
  ggplot(aes(x = as.factor(non.bloom.date), y = factor(bloom.date), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.96,0.98,1.00), 
                     labels=c("0.96","0.98","1.00"),
                     limits=c(.96,1.00))  +
  labs(y = "Bloom Dates", x = "All Dates") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust= 1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), face = "bold", size = 16),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"))

plot0608
plot0706
plot0720
plot0913
plot1026

plot2 <- compare.set.peaks %>%  
  ggplot(aes(x = as.factor(name1_date), y = factor(name2_date, levels = c("2012-10-26", "2012-09-13", "2012-07-20", "2012-07-06")), fill = conANI)) +
  geom_tile() +
  scale_fill_viridis(na.value = "gray",
                     breaks=c(.96,0.98,1.00), 
                     labels=c("0.96","0.98","1.00"),
                     limits=c(.96,1.00))  +
  labs(y = "Bloom Dates", x = "All Dates") +
  scale_x_discrete(labels=c('1', '2', '3', '4')) +
  scale_y_discrete(labels=c('5', '4', '3', '2')) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text=element_text(size=16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(1, 0.5, 1, 0.5), 
                           "inches"))

                                             