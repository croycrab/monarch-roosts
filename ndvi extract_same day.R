#download NDVI data from NCEI 

# load packages -----------------------------------------------------------
library(RCurl)
library(tidyverse)
library(stringr)
library(glue)
library(ncdf4)
library(sf)
library(terra)
library(exactextractr)

# load relevant files -----------------------------------------------------
load("GIS utilities.RData") 
map2$ST_12[c(73,63,64,65)] = "MX"
flyway = c("NLE","COA","TAM","TX","OK","MO","KS","AR","IA","IL","MI","MN","OH","IN","SD","ND","WI","ON","LA","NE")
map3 = map2 %>% filter(ST_12 %in% c(flyway,"MX"))

# get ncei_crs
ncei_crs = rast("ndvi raw/AVHRR-Land_v005_AVH13C1_NOAA-18_20070827_c20170401040337.nc") %>% st_crs()

#view circles
#view 12.5 km radius circles
read_csv("monarch roost_curated.csv") %>%
  filter(year >= 2007 & `State/Prov` %in% flyway) %>% 
  mutate(year = as.numeric(str_split(date,"-",simplify = T)[,1])) %>% 
  distinct(year,julian,lat,lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>%
  mutate(easting = st_coordinates(.)[, 1],
         northing = st_coordinates(.)[, 2],
         obs_id = as.integer(factor(paste0(lat,lon,year)))) %>%  
  st_buffer(dist = 12.5) %>% #25 km diameter
  ggplot()+
  geom_sf()+
  ggspatial::annotation_scale()

#get date by location unique combinations
sites = read_csv("monarch roost_curated.csv") %>%
  filter(year >= 2007 & `State/Prov` %in% flyway) %>% 
  mutate(year = as.numeric(str_split(date,"-",simplify = T)[,1])) %>% 
  distinct(year,julian,lat,lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>%
  mutate(easting = st_coordinates(.)[, 1],
         northing = st_coordinates(.)[, 2],
         obs_id = as.integer(factor(paste0(lat,lon,year)))) %>%  
  st_buffer(dist = 12.5) %>% #12.5 km radius
  st_transform(crs = ncei_crs)


# extract ndvi ------------------------------------------------------------
years = c(2007:2023)
ndvi_list_combined = list()
for(i in 1:length(years)){
  #get list of files for year
  ncei_url = glue("https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/{years[i]}/")
  site_filt = sites %>% 
    filter(year == years[i])
  if(years[i] < 2014){
    file_list <- getURL(ncei_url, ftp.use.epsv = FALSE, dirlistonly = TRUE) %>% 
      str_extract_all(">AVHRR.*\\.nc") %>% 
      unlist() %>% 
      unique() %>% 
      str_sub(start = 2) %>% 
      as.data.frame() %>% 
      slice(unique(site_filt$julian)) %>% 
      pull(1)
  } else {
    file_list <- getURL(ncei_url, ftp.use.epsv = FALSE, dirlistonly = TRUE) %>% 
      str_extract_all(">VIIRS.*\\.nc") %>% 
      unlist() %>% 
      unique() %>% 
      str_sub(start = 2) %>% 
      as.data.frame() %>% 
      slice(unique(site_filt$julian)) %>% 
      pull(1)
  }
  #download each file into temp directory
  site_filt = site_filt %>% 
    mutate(jul_id = as.integer(factor(julian)))
  ndvi_list = list()
  for(j in 1:length(unique(site_filt$jul_id))){
    site_filt_jul = site_filt %>% filter(jul_id == j) 
    url_temp = paste0(ncei_url,file_list[j])
    download.file(url = url_temp, destfile = paste0("ndvi raw/temp files/",file_list[j]))
    ndvi.rast = rast(paste0("ndvi raw/temp files/",file_list[1]))[["NDVI"]]
    ndvi.rast[ndvi.rast < .1] = NA
    ndvi_list[[j]] = data.frame(obs_id = site_filt_jul$obs_id,
                                julian = site_filt_jul$julian,
                                jul_id = site_filt_jul$jul_id,
                                year = years[i],
                                ndvi = exact_extract(ndvi.rast,site_filt_jul,'mean')
                                )
  }
  ndvi.df = do.call("rbind", ndvi_list)
  ndvi_list_combined[[i]] = ndvi.df
  unlink("ndvi raw/temp files/*")
}
ndvi.combined.df = do.call("rbind", ndvi_list_combined) %>% 
  left_join(sites %>% 
              st_drop_geometry() %>% 
              dplyr::select(obs_id,lat,lon),by = "obs_id")
write.csv(ndvi.combined.df,"ndvi_same day_25km.csv")

