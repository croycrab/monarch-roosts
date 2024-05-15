#curate data for analysis and inspect

#load libraries
library(openxlsx)
library(tidyverse)
library(sf)
library(ggpubr) #for pubclean theme

#import data, get maximum count per site*year
d1 = openxlsx::read.xlsx("monarch roost_original.xlsx") %>%
  rename(date = Date, lat = Lat, lon = Long,count = `#.monarchs`,year = Year, julian = Julian) %>% 
  mutate(date = as.character(as.Date(date,origin = "1899-12-30"))) %>% 
  filter(!(is.na(count)))  #remove NA observations

d1 %>% nrow() #4984 total observations
d1 %>% filter(year %in% 2003:2006) %>% nrow()

d1 = d1 %>% 
  group_by(lat,lon,year) %>% 
  mutate(count_max = max(count)) %>% 
  ungroup() %>% 
  distinct(lat,lon,year, .keep_all=TRUE) %>% 
  mutate(count = count_max) %>% 
  select(-count_max)
d1 %>% filter(year >= 2007) %>% nrow() #3080 total observations

write_csv(d1,"monarch roost_curated.csv")

#view number of observations over time
d1 %>% 
  group_by(year) %>% 
  summarize(n_obs = sum(count > 0)) %>% 
  ggplot(aes(y = n_obs,x = year)) +
  geom_point()+
  theme_pubclean() +
  labs(y = "Observations" , x = "Year")+
  geom_smooth(method = "loess",span = 0.5, se = F, color = "black")
ggsave("plots/appendix_observations by year.png",bg = "white",width = 5,height = 4)

#view distribution of observations over time
#get 2.5, 50, 97.5 quantiles
hist.quants = d1 %>% 
  mutate(log_count = log(count)) %>% 
  # filter(year >= 2007) %>% 
  pull(log_count) %>% 
  quantile(prob = c(0.025,0.5,0.975)) %>% 
  as.data.frame()

gghistogram(data = d1 %>%   
              # filter(year >= 2007) %>% 
              mutate(log_count = log(count)), 
            x = "log_count",y = "count",fill= "gray")+
  geom_vline(xintercept = hist.quants[1,1],lty = "dashed")+ #2.5%
  geom_vline(xintercept = hist.quants[2,1],lty = "dashed")+ #median
  geom_vline(xintercept = hist.quants[3,1],lty = "dashed")+ #97.5%
  facet_wrap(~year)+
  theme_pubclean()+
  labs(y = "Observations", x = "ln-Roost Size")
ggsave("plots/appendix_histogram of observations by year_all years.png",bg = "white",width = 10,height = 8)

#create another cleaned up version for potential inclusion as main figure
gghistogram(data = d1 %>%   
              filter(year %in% c(2008,2011,2014,2017,2020,2023)) %>%
              mutate(log_count = log(count)), 
            x = "log_count",y = "count",fill= "gray")+
  geom_vline(xintercept = hist.quants[1,1],lty = "dashed")+ #2.5%
  geom_vline(xintercept = hist.quants[2,1],lty = "dashed")+ #median
  geom_vline(xintercept = hist.quants[3,1],lty = "dashed")+ #97.5%
  facet_wrap(~year)+
  theme_pubclean()+
  labs(y = "Observations", x = "ln-Roost Size")
ggsave("plots/histogram of observations by year_6 panels.png",bg = "white",width = 10,height = 6)


#view map of distributions of observations over time
## create and save some GIS utilities 
#GIS Utilities -- create map and define CRS
#load
#define a CRS -- using CRS from Meehan et al (2023 -- SVCs using inlabru tutorial)
#this CRS should work fine, 
epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

#create sites sf object
sites <- d1 %>% 
  distinct(lat, lon) %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs = epsg6703km)

#make a base map
#I like how easy it is to download the US/CA shapefile from the bbsBayes package, so I'll do that
map1 <- bbsBayes::load_map(stratify_by="state") %>% 
  st_transform(crs=epsg6703km) #reproject vector to common CRS
#need to add Mexico
mexico = st_read("mexico shapefile/mexican-states.shp") %>% 
  st_transform(crs=st_crs(map1)) #reproject vector to common CRS

#add mexico to us/ca map, cropping with a buffer for well-sized plots later
map2 = map1 %>% 
  bind_rows(mexico %>% select(geometry)) %>% 
  st_crop(st_bbox(sites %>% st_buffer(dist = 500)))

save(epsg6703km,map2,file = "GIS utilities.RData")

d1 %>% 
  distinct(lat, lon, year) %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs = epsg6703km) %>%
  ggplot() +
  theme_bw()+
  geom_sf(inherit.aes = F, data = map2 %>% 
            st_simplify(preserveTopology = F,dTolerance = 10)) +
  geom_sf(alpha=0.8,size = 1.25,pch=4)+
  facet_wrap(~year)
ggsave("plots/appendix_observation map by year.png",bg = "white",height = 8,width = 12)



