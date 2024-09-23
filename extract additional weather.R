#scrape weather data from additional sources (DAYMET, CRU-TS) and compare with NOAA station weather data. This is in response to reviewer feedback that more localized weather phenomenon might be overlooked if we are using weighted averages of weather station data within a 50km radius of an observation. There are different interpolation methods out there that have been developed to fill this gap, and so I'll grab data from those sources and compare with our "cruder" approach to see if they're more or less congruent. 

#import observation info (time and place)
load("weather.rda")
df.week2 = df.week %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326)
#PRISM
library(prism)
library(sf)
library(terra)
library(stringr)
library(exactextractr)
library(tidyverse)

unique_dates = unique(df.week$date)
clim_vars = c("ppt","tmean")
prism_set_dl_dir("PRISM/")

ppt.list = list()
tmean.list = list()
for(i in 1:length(unique_dates)){
  for(j in 1:length(clim_vars)){
    get_prism_dailys(type = clim_vars[j],dates = unique_dates[i])
    file_name = paste0("PRISM_",clim_vars[j], "_stable_4kmD2_",str_replace_all(unique_dates[i], "-",""),"_bil")
    assign(paste0(clim_vars[j],".rast"),rast(paste0("PRISM/",file_name,"/",file_name,".bil")))
  }
  ppt.list[[i]] =  extract(ppt.rast,df.week2 %>% 
            filter(date == unique_dates[i]) %>% 
            st_transform(., crs = st_crs(ppt.rast))) %>% 
    rename(ppt = 2) %>% 
    mutate(ID = df.week2 %>% 
             filter(date == unique_dates[i]) %>% 
             pull(id),
           date = unique_dates[i])
  
  tmean.list[[i]] =  extract(tmean.rast,df.week2 %>% 
                             filter(date == unique_dates[i]) %>% 
                             st_transform(., crs = st_crs(tmean.rast))) %>% 
    rename(tmean = 2) %>% 
    mutate(ID = df.week2 %>% 
             filter(date == unique_dates[i]) %>% 
             pull(id),
           date = unique_dates[i])

}
ppt.df = do.call("rbind",ppt.list)
tmean.df = do.call("rbind",tmean.list)

prism.df = ppt.df %>% 
  left_join(tmean.df, by = c("ID","date")) %>% 
  left_join(df.week,  by = c("ID" = "id","date"))

write.csv(prism.df,"prism.csv")
# unlink("PRISM/*")

#DAYMET
library(daymetr)
daymet.list = list()
daymet.ids = unique(df.week$id)
for(i in 1:length(daymet.ids)){
  lat_i = df.week %>% 
    filter(id == i) %>% 
    pull(lat) %>% 
    unique()
  lon_i = df.week %>% 
    filter(id == i) %>% 
    pull(lon) %>% 
    unique()
  year_i = df.week %>% 
    filter(id == i) %>% 
    pull(date) %>% 
    str_sub(.,1,4) %>% 
    unique() %>% 
    as.numeric()
  daymet.list[[i]] = try(
    download_daymet(
      lat = lat_i,
      lon = lon_i,
      start = year_i,
      end = year_i
      )
  )
}
daymet.df = data.frame(year = NA, yday = NA, precip = NA, tmax = NA, tmin = NA)
for(i in 1:length(daymet.ids)){
  if((daymet.list[[i]][1] %>% grep(.,pattern = "Error") %>% is_empty())){
    daymet.df = bind_rows(daymet.df,daymet.list[[i]]$data %>% 
                            rename(precip = prcp..mm.day.,
                                   tmax = tmax..deg.c.,
                                   tmin = tmin..deg.c.) %>% 
                            select(year, yday, precip, tmax, tmin) %>% 
                            mutate(id = daymet.ids[i],
                                   date = as.Date(paste(year, yday, sep = "-"), "%Y-%j")) %>% 
                            filter(date %in% (df.week %>% filter(id == daymet.ids[i]) %>% pull(date)))
    )
  }
    # else(
    #      )
}
write.csv(daymet.df %>% na.omit(),"daymet.csv")


clim_vars2 = c("pre","tmp")
tmp.rast = rast("CRU-TS/cru_ts4.08.2001.2010.tmp.dat.nc")
tmp.rast

#drought index
#from: https://droughtmonitor.unl.edu/DmData/DataDownload/ComprehensiveStatistics.aspx
drought.df = read_csv("drought index_weekly.csv") %>% 
  mutate(year = str_sub(MapDate, 1,4) %>% as.numeric(),
         month = str_sub(MapDate, 5,6) %>% as.numeric()) %>% 
  filter(month %in% 8:11) %>% 
  pivot_longer(cols = None:D4,names_to = "drought", values_to = "area") %>% 
  group_by(year, drought,Region) %>% 
  reframe(area = mean(area)) %>% 
  ungroup()


ggplot(drought.df %>% mutate(drought = fct_relevel(drought,c("None","D0","D1","D2","D3","D4"))), aes(y = area, x = year,group = Region, color = Region))+
  geom_point()+
  geom_smooth(method = "lm", se = T)+
  ggpubr::theme_pubclean()+
  facet_grid(~drought)+
  labs(y = "Area Affected by Drought", x = "Year", color = "Region")
ggsave("plots/area affected by drought.png", width = 10, height = 4)

mod.none = lm(area~year*Region, data = drought.df %>% filter(drought == "None"))
car::Anova(mod.none, type = 2)

mod.d0 = lm(area~year*Region, data = drought.df %>% filter(drought == "D0"))
car::Anova(mod.d0, type = 2)

mod.d1 = lm(area~year*Region, data = drought.df %>% filter(drought == "D1"))
car::Anova(mod.d1, type = 2)

mod.d2 = lm(area~year*Region, data = drought.df %>% filter(drought == "D2"))
car::Anova(mod.d2, type = 2)

mod.d3 = lm(area~year*Region, data = drought.df %>% filter(drought == "D3"))
car::Anova(mod.d3, type = 2)

mod.d4 = lm(area~year*Region, data = drought.df %>% filter(drought == "D4"))
car::Anova(mod.d4, type = 2)
