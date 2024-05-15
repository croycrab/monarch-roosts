#prepare weather csv for analysis

#load libraries
library(tidyverse)

load("weather.rda")
vcweather.df = read.csv("vcweather.csv") %>% 
  mutate(date = as.Date(date)) %>% 
  select(-X)

df.week2 = df.week %>% 
  left_join(vcweather.df,by = c("date","lat","lon"), relationship = "many-to-many") %>% 
  distinct(lat, lon, date, id, .keep_all = T) %>% 
  # filter(lat > 26) %>% 
  group_by(lat, lon, id) %>% 
  mutate(winddir_full = winddir,
         winddir = abs(winddir - 180),
         moonphase = abs(moonphase - 0.5)) %>% 
  summarise(pressure_24hr = pressure[which(date == max(date))]-pressure[which(date == (max(date)-1))],
            temp7 = mean(temp, na.rm = T),
            precip7 = sum(precip,na.rm = T),
            winddir7 = mean(winddir, na.rm = T),
            windspeed7 = mean(windspeed, na.rm = T),
            winddir7_full = mean(winddir_full, na.rm = T),
            date = max(date))

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
df.week2[is.nan(df.week2)] <- NA

is.inf.data.frame <- function(x)
  do.call(cbind, lapply(x, is.infinite))
df.week2[is.inf.data.frame(df.week2)] <- NA
  
# df.week2 = df.week2[complete.cases(df.week2),]  
write.csv(df.week2,"curated weather.csv")



