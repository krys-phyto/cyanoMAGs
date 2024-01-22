# Krystyn Kibler 2023-01-11
# Purpose: Figure that shows temp with depth and sample date lines


### Libraries ###
library(rio)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(data.table)
library(ggpubr)
library(egg)


### Files ###
environment <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/environmental/Mendota-envdata-2008-2012-edit.xlsx"

file_relabund <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/anvio/anvi-metagenomics/Refined-SUMMARY/bins_across_samples/abundance.txt"

# reading data from all sheets
data <- import_list(environment)

# print data
print (data)

# need only temp data 
temp.buoy <- data$`WaterTemp-buoy`
temp.mo <- data$`WaterTemp-YSI`

# relabund
relabund <- read_delim(file_relabund,
                       "\t", escape_double = FALSE, trim_ws = TRUE)


### Temperature Figure ###

### Metag Sampledates ###
  # tidy data
  # current header format = metag_20080719_OR WANT just the date
metag_header_dates <- colnames(relabund)
metag_header_dates <- sapply(strsplit(metag_header_dates, "_"), "[", 2)
metag_header_dates <- as_date(metag_header_dates) # gets metagdates too

sampledates <- as.data.frame(metag_header_dates)

### MO temp ###
  #tidy data
  # separate sample date and time into two columns
temp.mo <- temp.mo %>% separate(Date, c("sampledate","time"), sep = " ")

  # remove unneccesary columns and change character na's into true na's
temp.mo <- temp.mo %>% select(-c(26:31))
temp.mo[ temp.mo == "NA" ] <- NA

  # turn all columns with values into numeric
temp.mo <- temp.mo %>% 
  mutate_at(c(3:25), as.numeric)

temp.mo <- temp.mo %>% pivot_longer(c(3:25), values_to = "wtemp", names_to = "depth")

  # convert sampledate column to date with lubridate
temp.mo$sampledate <- as_date(temp.mo$sampledate)

  # convert depth column into numeric
temp.mo$depth <- as.numeric(temp.mo$depth)

  # remove rows where there is an na in the sampledate column
temp.mo <- temp.mo %>% 
  drop_na(c("sampledate"))

### Data Interpolations ###
  # data interpolations through dates #
estimate_temp_by_date <- function(target_date, target_depth) {
  data_for_date <- temp.mo %>% 
    filter(sampledate == target_date) %>%
    arrange(depth)
  
  # approx() is one way to do a linear interpolation
  approx(data_for_date$depth, data_for_date$wtemp, xout = target_depth)$y
}

estimate_temp_by_date(ymd("2008-06-02"), c(5.0, 6.0, 7.0, 8.0))
estimate_temp_by_date(ymd("2008-06-02"), c(5.5, 6.5, 7.5))

temp_interp_depth <- crossing(
  # the same dates as sonde_tbl_1993
  tibble(date = unique(temp.mo$sampledate)),
  # depths can now be any value
  tibble(depth = seq(0, 20, length.out = 150))
) %>%
  group_by(date) %>%
  mutate(temp = estimate_temp_by_date(date[1], depth))


  # data interpolations through time #
estimate_temp_by_depth <- function(target_depth, target_date) {
  data_for_depth <- temp_interp_depth %>% 
    filter(depth == target_depth) %>%
    arrange(date)
  approx(data_for_depth$date, data_for_depth$temp, xout = target_date)$y
}

estimate_temp_by_depth(
  target_depth = 0, 
  target_date = seq(ymd("2008-06-02"), ymd("2008-06-09"), by = 1)
)

  # temp interpolations
  # 2008 #
temp_raster_2008 <- crossing(
  # dates can now be any value
  tibble(date = seq(ymd("2008-06-02"), ymd("2008-10-23"), by = 1)),
  # depths must be the same as in temp_interp_depth
  tibble(depth = unique(temp_interp_depth$depth))
) %>%
  group_by(depth) %>%
  mutate(temp = estimate_temp_by_depth(depth[1], date))

  # 2009 #
temp_raster_2009 <- crossing(
  tibble(date = seq(ymd("2009-04-10"), ymd("2009-10-17"), by = 1)),
  tibble(depth = unique(temp_interp_depth$depth))
) %>%
  group_by(depth) %>%
  mutate(temp = estimate_temp_by_depth(depth[1], date))

  # 2010 #
temp_raster_2010 <- crossing(
  tibble(date = seq(ymd("2010-05-19"), ymd("2010-08-31"), by = 1)),
  tibble(depth = unique(temp_interp_depth$depth))
) %>%
  group_by(depth) %>%
  mutate(temp = estimate_temp_by_depth(depth[1], date))

  # 2011 #
temp_raster_2011 <- crossing(
  tibble(date = seq(ymd("2011-05-06"), ymd("2011-09-17"), by = 1)),
  tibble(depth = unique(temp_interp_depth$depth))
) %>%
  group_by(depth) %>%
  mutate(temp = estimate_temp_by_depth(depth[1], date))

  # 2012 #
temp_raster_2012 <- crossing(
  tibble(date = seq(ymd("2012-06-08"), ymd("2012-11-09"), by = 1)),
  tibble(depth = unique(temp_interp_depth$depth))
) %>%
  group_by(depth) %>%
  mutate(temp = estimate_temp_by_depth(depth[1], date))

### Sampledate prep for plot ###
  # tidy
sampledates1 <- sampledates
sampledates2 <- sampledates

sampledates1['depth']=0
sampledates2['depth']=12

sampledates1 <- sampledates1[-1,]
sampledates2 <- sampledates2[-1,]

datlist <- list(sampledates1, sampledates2)
sampledates <- rbindlist(datlist)

sampledates$value <- 1

colnames(sampledates) <- c("date","depth", "value")


### Plotting even ###

# Modify temprasters to start at same month and date #
# Latest date in time series 2012-11-09
# Earliest date in time series 2009-04-30

# 2008 
temp_raster_2008[nrow(temp_raster_2008) + 1,] = list(date=as_date("2008-04-15"), depth.x=0, value=100/0)
temp_raster_2008[nrow(temp_raster_2008) + 1,] = list(date=as_date("2008-11-15"), depth.x=0, value=100/0)
temp_raster_2008[sapply(temp_raster_2008, is.infinite)] <- NA

# 2009
temp_raster_2009[nrow(temp_raster_2009) + 1,] = list(date=as_date("2009-04-15"), depth.x=0, value=100/0)
temp_raster_2009[nrow(temp_raster_2009) + 1,] = list(date=as_date("2009-11-15"), depth.x=0, value=100/0)
temp_raster_2009[sapply(temp_raster_2009, is.infinite)] <- NA

# 2010
temp_raster_2010[nrow(temp_raster_2010) + 1,] = list(date=as_date("2010-04-15"), depth.x=0, value=100/0)
temp_raster_2010[nrow(temp_raster_2010) + 1,] = list(date=as_date("2010-11-15"), depth.x=0, value=100/0)
temp_raster_2010[sapply(temp_raster_2010, is.infinite)] <- NA

# 2011
temp_raster_2011[nrow(temp_raster_2011) + 1,] = list(date=as_date("2011-04-15"), depth.x=0, value=100/0)
temp_raster_2011[nrow(temp_raster_2011) + 1,] = list(date=as_date("2011-11-15"), depth.x=0, value=100/0)
temp_raster_2011[sapply(temp_raster_2011, is.infinite)] <- NA

# 2012
temp_raster_2012[nrow(temp_raster_2012) + 1,] = list(date=as_date("2012-04-15"), depth.x=0, value=100/0)
temp_raster_2012[nrow(temp_raster_2012) + 1,] = list(date=as_date("2012-11-15"), depth.x=0, value=100/0)
temp_raster_2012[sapply(temp_raster_2012, is.infinite)] <- NA


temp_raster_2008 <- left_join(temp_raster_2008, sampledates, by = "date")
temp_raster_2009 <- left_join(temp_raster_2009, sampledates, by = "date")
temp_raster_2010 <- left_join(temp_raster_2010, sampledates, by = "date")
temp_raster_2011 <- left_join(temp_raster_2011, sampledates, by = "date")
temp_raster_2012 <- left_join(temp_raster_2012, sampledates, by = "date")






# 2008 plotting
plot_2008 <- temp_raster_2008 %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,25,15), name = "\n Temp \n (\u00B0C)") + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  labs(y = " ", title = " ") +
  guides(color = "none") +
  scale_x_date(date_labels="%b", date_breaks  ="1 month", breaks = seq.Date(from = as.Date("2008-04-15"), 
                                                                            to = as.Date("2008-11-15"), 
                                                                            by = "1 months")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(), 
        legend.text = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 14))


# 2009 plotting
plot_2009 <- temp_raster_2009 %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  scale_x_date(date_labels="%b", date_breaks  ="1 month", breaks = seq.Date(from = as.Date("2009-04-15"), 
                                                                            to = as.Date("2009-11-15"), 
                                                                            by = "1 months")) +
  labs(y = " ", title = "2008") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(), 
        title = element_blank(),
        legend.text = element_text(face = "bold", size = 14))

# 2010 plotting
plot_2010 <- temp_raster_2010 %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  scale_x_date(date_labels="%b", date_breaks  ="1 month", breaks = seq.Date(from = as.Date("2010-04-15"), 
                                                                            to = as.Date("2010-11-15"), 
                                                                            by = "1 months")) +
  ylab("Depth (m)") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(), 
        title = element_blank(),
        legend.text = element_text(face = "bold", size = 14))

# 2011 plotting
plot_2011 <- temp_raster_2011 %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  scale_x_date(date_labels="%b", date_breaks  ="1 month", breaks = seq.Date(from = as.Date("2011-04-15"), 
                                                                            to = as.Date("2011-11-15"), 
                                                                            by = "1 months")) +
  labs(y = " ", title = "2008") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        title = element_blank(),
        legend.text = element_text(face = "bold", size = 14))

# 2012 plotting
plot_2012 <- temp_raster_2012 %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  scale_x_date(date_labels="%b", date_breaks  ="1 month", breaks = seq.Date(from = as.Date("2012-04-15"), 
                                                                            to = as.Date("2012-11-15"), 
                                                                            by = "1 months")) +
  labs(y = " ", x = "Month", title = "2008") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_text(face = "bold", size = 14), 
        title = element_blank(),
        legend.text = element_text(face = "bold", size = 14))

plot_even <- ggarrange(plot_2008, plot_2009, plot_2010, plot_2011, plot_2012, 
                  ncol = 1, common.legend = TRUE, legend = "right", align = "h",
                  labels = c("2008", "2009", "2010", "2011", "2012"))










### Plotting uneven ###
# 2008 plotting
temp_raster_2008_uneven <- left_join(temp_raster_2008, sampledates, by = "date")

plot_2008 <- temp_raster_2008_uneven %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,25,15), name = "\n Temp \n (\u00B0C)") + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  labs(y = " ", title = "2008") +
  guides(color = "none") +
  theme(legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        axis.text.x=element_text(size=10), 
        axis.title.x=element_blank())

leg <- get_legend(plot_2008)
temp.leg <- as_ggplot(leg)

# 2009 plotting
temp_raster_2009_uneven <- left_join(temp_raster_2009, sampledates, by = "date")

plot_2009 <- temp_raster_2009_uneven %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size = 1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  theme(legend.position="none", 
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10)) +
  labs(x ="Date", y = " ", title = "2009")

# 2010 plotting
temp_raster_2010_uneven <- left_join(temp_raster_2010, sampledates, by = "date")

plot_2010 <- temp_raster_2010_uneven %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  theme(legend.position="none", 
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10)) +
  labs(x ="Date", y = "Depth", title = "2010")

# 2011 plotting
temp_raster_2011_uneven <- left_join(temp_raster_2011, sampledates, by = "date")

plot_2011 <- temp_raster_2011_uneven %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  theme(legend.position="none", 
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10)) +
  labs(x ="Date", y = " ", title = "2011")

# 2012 plotting
temp_raster_2012_uneven <- left_join(temp_raster_2012, sampledates, by = "date")

plot_2012 <- temp_raster_2012_uneven %>% ggplot() +
  geom_raster(aes(date, depth.x, fill = temp)) +
  scale_y_reverse() +
  scale_fill_viridis(option = "D", breaks = c(10,20,10)) + 
  coord_cartesian(expand = FALSE) +
  geom_line(aes(date, depth.y, color = value, group = date), size=1.25) +
  scale_colour_gradient(low = "black", high = "black") +
  theme(legend.position="none",
        axis.text.x=element_text(size=10)) +
  labs(x ="Date", y = " ", title = "2012")

plot_uneven <- ggarrange(plot_2008, plot_2009, plot_2010, plot_2011, plot_2012,
                         nrow = 5, ncol = 1)



