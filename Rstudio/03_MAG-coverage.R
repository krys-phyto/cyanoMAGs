# Krys Kibler
# Purpose: Coverage between MAGs


### Files ###
file_relabund <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bins_across_samples/abundance.txt"
file_nicknames <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/mags/bins_nicknames.txt"
file_taxonomy <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvio-taxonomy/SCG-TAXONOMY.txt"


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
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
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


### Environmental parameters surrounding MAG detection days 




