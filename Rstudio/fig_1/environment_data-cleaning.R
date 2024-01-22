# Krys Kibler 
# 2023-10-06
# Environmental Data


### Libraries ###
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(stringr)
library(dplyr)
library(scales)

### Data ###
DO <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/environmental/DO.xlsx"
nitrate <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/environmental/nitrate.xlsx"
pH <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/environmental/pH.xlsx"
SRP <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/environmental/SRP.xlsx"
temp <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/paper-1-cyanoMags/data/environmental/ysi_temp.xlsx"

DO <- read_excel(DO)
nitrate <- read_excel(nitrate)
pH <- read_excel(pH)
SRP <- read_excel(SRP)
temp <- read_excel(temp)

### Tidy ###

# DO
DO <- DO %>% select(c(1:19))
DO$Date <- gsub(" ..:..:..", "", DO$Date)
DO$Date <- as_date(DO$Date)
DO <- DO %>% mutate_at(c(2:19), as.numeric)
DO <- DO %>% pivot_longer(c(2:19), names_to = "depth", values_to = "DO")
DO$year4 <- year(DO$Date)


# nitrate
nitrate <- nitrate %>% select(-3)
nitrate$Date <- gsub(" ..:..:..", "", nitrate$Date)
nitrate$Date <- as_date(nitrate$Date)
nitrate$year4 <- year(nitrate$Date)

# pH
pH <- pH %>% select(-3)
pH$sampledate <- as_date(pH$sampledate)
pH$year4 <- year(pH$sampledate)

# SRP 
SRP <- SRP %>% select(-3)
SRP$Date <- gsub(" ..:..:..", "", SRP$Date)
SRP$Date <- as_date(SRP$Date)
SRP$year4 <- year(SRP$Date)

# temp
temp.12 <- temp %>% select(c(1:16))
temp.12$Date <- gsub(" ..:..:..", "", temp.12$Date)
temp.12$Date <- as_date(temp.12$Date)
temp.12 <- temp.12 %>% mutate_at(c(2:16), as.numeric)
temp.12 <- temp.12 %>% pivot_longer(c(2:16), names_to = "depth", values_to = "temp")
temp.12$year4 <- year(temp.12$Date)


### calculations ###
# DO
DO <- DO %>% filter(DO < 200)

DO.year <- DO %>% group_by(year4) %>% 
  summarise(mean=mean(DO), sd=sd(DO))

# nitrate
nitrate.year <- nitrate %>% group_by(year4) %>% 
  summarise(mean=mean(`NO3 (mgN/L)`), sd=sd(`NO3 (mgN/L)`))

# pH
pH.year <- pH %>% group_by(year4) %>% 
  summarise(mean=mean(pH), sd=sd(pH))

# SRP
SRP.year <- SRP %>% group_by(year4) %>% 
  summarise(mean=mean(`DRP (mg-P/L)`), sd=sd(`DRP (mg-P/L)`))

# temp.12
temp.12.year <- temp.12 %>% drop_na() %>% group_by(year4) %>% 
  summarise(mean=mean(temp), sd=sd(temp))



### Figure 1

temp <- temp %>% select(c(1:24))
temp$Date <- gsub(" ..:..:..", "", temp$Date)
temp$Date <- as_date(temp$Date)
temp <- temp %>% mutate_at(c(2:24), as.numeric)
temp <- temp %>% pivot_longer(c(2:24), names_to = "depth", values_to = "temp")
temp$year4 <- year(temp$Date)

fig1 <- temp %>% drop_na() %>% filter(year4 == 2008) %>% 
  ggplot(aes(x = Date, y = depth)) +
  geom_raster(aes(fill = temp), interpolate = TRUE)


geom_raster(aes(fill = Temp, x = Longitude), interpolate = TRUE) +
  




