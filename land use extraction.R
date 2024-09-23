#extract urban area
#land use
#load packages
library(dplyr); library(raster); library(sf); library(bbsBayes); library(ncdf4); library(tidyr); library(stringr); library(RCurl); library(XML); library(ncdf4); library(terra); library(exactextractr)

#set wd
#load roost observations and GIS utilities
load("roost data.rda")
load("GIS utilities.RData") 

#import land use data
# lulc.rast <- rast("~/Dropbox/NABA2/land use/hildap_vGLOB-1.0-f_netcdf/hildaplus_GLOB-1-0-f_states.nc")
lulc.rast <- rast("land use data/hildap_vGLOB-1.0-f_netcdf/hildaplus_GLOB-1-0-f_states.nc")

#view 5 km radius circles
count_dat %>% 
  distinct(lat, lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>% 
  st_buffer(dist = 5) %>% 
  ggplot()+
  geom_sf()+
  ggspatial::annotation_scale()
  
#create 5km radius circles around observations and clip raster
circ.df = count_dat %>% 
  distinct(lat, lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>% 
  st_buffer(dist = 10) %>% 
  st_transform(st_crs(lulc.rast))

lulc.rast2 = lulc.rast[[109:121]] %>% 
  crop(.,ext(circ.df)) %>% 
  app(., fun=function(x){ x[x < 11] <- 0; return(x)}) %>% 
  app(., fun=function(x){ x[x > 11] <- 0; return(x)}) %>% 
  classify(., cbind(id = c(11), v = c(1)))

lulc.list = list()
yr.vec = seq(2007,2019,1)
for(i in 1:nlyr(lulc.rast2)){
  lulc.list[[i]] = exact_extract(x=lulc.rast2[[i]],y=circ.df,fun = "mean") %>% 
    as.data.frame() %>% 
    rename(urban = 1) %>% 
    mutate(circ_id = 1:nrow(circ.df),
           year = yr.vec[i],
           ) %>% 
    dplyr::select(circ_id,year,urban)
}
lulc.df <- do.call("rbind",lulc.list)

lulc.df2 = lulc.df %>% 
  group_by(circ_id) %>% 
  summarize(urban = mean(urban))

png("plots/urban habitat around observations.png")
hist(lulc.df2$urban, breaks = 25,xlab = "Proportion of Urban Habitat",main = NULL)
dev.off()

lulc.df2 %>% 
  filter(urban < 0.1) %>% 
  pull(urban) %>% 
  length()
1425/nrow(lulc.df2)
write.csv(lulc.df2,"urban habitat.csv")

#extract natural habitat
m = c(0,39,0,
      56,77,0,
      40,55,1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
lulc.rast3 = lulc.rast[[109:121]] %>% 
  crop(.,ext(circ.df)) %>% 
  classify(., rclmat)

lulc.list3 = list()
yr.vec = seq(2007,2019,1)
for(i in 1:nlyr(lulc.rast3)){
  lulc.list3[[i]] = exact_extract(x=lulc.rast3[[i]],y=circ.df,fun = "mean") %>% 
    as.data.frame() %>% 
    rename(natural = 1) %>% 
    mutate(circ_id = 1:nrow(circ.df),
           year = yr.vec[i],
    ) %>% 
    dplyr::select(circ_id,year,natural)
}
lulc.df3 <- do.call("rbind",lulc.list3)
lulc.df4 = lulc.df3 %>% 
  group_by(circ_id) %>% 
  summarize(natural = mean(natural))

png("plots/natural habitat around observations.png")
hist(lulc.df4$natural, breaks = 25,xlab = "Proportion of Natural Habitat",main = NULL)
dev.off()

#view whether the sites with at least two observations (which generate the trends) differ in urban habitat
circ.df_pair = count_dat %>% 
  distinct(lat, lon, year) %>%
  group_by(lat,lon) %>% 
  mutate(n_obs = n()) %>% 
  ungroup() %>% 
  filter(n_obs >= 2) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>% 
  st_buffer(dist = 10) %>% 
  st_transform(st_crs(lulc.rast))

lulc.list2 = list()
yr.vec = seq(2007,2019,1)
for(i in 1:nlyr(lulc.rast2)){
  lulc.list2[[i]] = exact_extract(x=lulc.rast2[[i]],y=circ.df_pair,fun = "mean") %>% 
    as.data.frame() %>% 
    rename(urban = 1) %>% 
    mutate(circ_id = 1:nrow(circ.df_pair),
           year = yr.vec[i],
    ) %>% 
    dplyr::select(circ_id,year,urban)
}
lulc.df3 <- do.call("rbind",lulc.list2)

lulc.df4 = lulc.df3 %>% 
  group_by(circ_id) %>% 
  summarize(urban = mean(urban))

hist(lulc.df4$urban, breaks = 25,xlab = "Proportion of Urban Habitat",main = NULL)

lulc.df4 %>% 
  filter(urban < 0.1) %>% 
  pull(urban) %>% 
  length()
794/nrow(lulc.df4)
