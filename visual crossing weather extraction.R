#download data from visual crossing
api_key = NA

# daily weather extraction from weather stations

#load libraries
library(tidyverse); library(sf); library(glue); library(lubridate); library(rnoaa);
library(rjson)

#create lat/lon with dates up to 7 days prior to an observation for extraction
lat_lon.df = read_csv("monarch roost_curated.csv") %>%
  distinct(date,lat,lon,year) %>% 
  mutate(id = row_number())

id.list = list()
for(i in 1:nrow(lat_lon.df)){
  id.df = lat_lon.df %>% 
    filter(id == i)
  id.df2 = id.df
  for(j in 1:6){
    id.df2 = id.df2 %>% 
      rbind(.,id.df %>% 
              mutate(date = date - j))
  }
  id.list[[i]] = id.df2
}
df.week =  do.call("rbind", id.list) %>% 
  arrange(id)
save(df.week,file = "weather.rda")

#get weather data
vcweather.list = list()
for(i in 1:nrow(df.week)){
  vc_url = glue('https://weather.visualcrossing.com/VisualCrossingWebServices/rest/services/timeline/{paste0(df.week$lat[i],",",df.week$lon[i])}/{df.week$date[i]}?key={api_key}')
  json_data <- fromJSON(file=vc_url)
  vcweather.list[[i]] = json_data$days %>% 
    do.call("rbind",.) %>%  
    as.data.frame() %>%  
    mutate(lat = df.week$lat[i], lon = df.week$lon[i]) %>% 
    rename(date = datetime)
  closeAllConnections()
}
vcweather.df = bind_rows(vcweather.list)
save(vcweather.df,file = "vcweather.rda")

#elements within the df are lists (?), so im going to run through with a for loop to extract numerics...
load("vcweather.rda")
vcweather.df2 = vcweather.df %>% 
  select(lat,lon)
for(i in 1:nrow(vcweather.df2)){
  vcweather.df2$date[i] = ifelse(is.null(vcweather.df$date[[i]]),NA,vcweather.df$date[[i]])
  vcweather.df2$tempmin[i] = ifelse(is.null(vcweather.df$tempmin[[i]]),NA,vcweather.df$tempmin[[i]])
  vcweather.df2$tempmax[i] = ifelse(is.null(vcweather.df$tempmax[[i]]),NA,vcweather.df$tempmax[[i]])
  vcweather.df2$temp[i] = ifelse(is.null(vcweather.df$temp[[i]]),NA,vcweather.df$temp[[i]])
  vcweather.df2$feelslikemax[i] = ifelse(is.null(vcweather.df$feelslikemax[[i]]),NA,vcweather.df$feelslikemax[[i]])
  vcweather.df2$feelslikemin[i] = ifelse(is.null(vcweather.df$feelslikemin[[i]]),NA,vcweather.df$feelslikemin[[i]])
  vcweather.df2$feelslike[i] = ifelse(is.null(vcweather.df$feelslike[[i]]),NA,vcweather.df$feelslike[[i]])
  vcweather.df2$dew[i] = ifelse(is.null(vcweather.df$dew[[i]]),NA,vcweather.df$dew[[i]])
  vcweather.df2$humidity[i] = ifelse(is.null(vcweather.df$humidity[[i]]),NA,vcweather.df$humidity[[i]])
  vcweather.df2$precip[i] = ifelse(is.null(vcweather.df$precip[[i]]),NA,vcweather.df$precip[[i]])
  vcweather.df2$windgust[i] = ifelse(is.null(vcweather.df$windgust[[i]]),NA,vcweather.df$windgust[[i]])
  vcweather.df2$windspeed[i] = ifelse(is.null(vcweather.df$windspeed[[i]]),NA,vcweather.df$windspeed[[i]])
  vcweather.df2$winddir[i] = ifelse(is.null(vcweather.df$winddir[[i]]),NA,vcweather.df$winddir[[i]])
  vcweather.df2$pressure[i] = ifelse(is.null(vcweather.df$pressure[[i]]),NA,vcweather.df$pressure[[i]])
  vcweather.df2$cloudcover[i] = ifelse(is.null(vcweather.df$cloudcover[[i]]),NA,vcweather.df$cloudcover[[i]])
  vcweather.df2$moonphase[i] = ifelse(is.null(vcweather.df$moonphase[[i]]),NA,vcweather.df$moonphase[[i]])
  vcweather.df2$source[i] = ifelse(is.null(vcweather.df$source[[i]]),NA,vcweather.df$source[[i]])
}
str(vcweather.df2)
write_csv(vcweather.df2,"vcweather.csv")
