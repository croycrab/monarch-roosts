#build models per the approach outlined in Meehan et al (2023 -- SVCs using SPDE in inlabru)

# load libraries ----------------------------------------------------------
#load libraries
library(maps)
library(tidyverse)
library(sf)
library(terra)
library(tidyterra) # raster plotting
library(scales)
library(INLA)
library(inlabru)
library(fmesher)
library(ggpubr)
library(ggspatial)

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
# Note: the 'splancs' package also needs to be installed,
# but doesn't need to be loaded

# set option
select <- dplyr::select
options(scipen = 99999)
options(max.print = 99999)
options(stringsAsFactors = FALSE)


# import data -------------------------------------------------------------
#maps / CRS
load("GIS utilities.RData") 

#adjust map
map2$ST_12[c(73,63,64,65)] = "MX"
flyway = c("NLE","COA","TAM","TX","OK","MO","KS","AR","IA","IL","MI","MN","OH","IN","SD","ND","WI","ON","LA","NE")
map3 = map2 %>% filter(ST_12 %in% c(flyway,"MX"))

#import data
weather.df = read_csv("curated weather.csv") %>% 
  select(date,lat,lon,temp7,precip7,winddir7,windspeed7,winddir7_full) %>%
  mutate(date = as.character(date),
         year = as.numeric(str_split(date,"-",simplify = T)[,1]),
         temp = temp7,
         rain = precip7,
         tailwind = windspeed7*cos(winddir7_full*(pi/180)))
ndvi.df = read_csv("ndvi_same day.csv")


count_dat = read_csv("monarch roost_curated.csv") %>%
  # select(-1) %>% #remove first column 
  filter(year >= 2007 & `State/Prov` %in% flyway) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>%
  mutate(easting = st_coordinates(.)[, 1],
         northing = st_coordinates(.)[, 2]) %>% 
  # dplyr::right_join(weather.df, by = c("date","lat","lon")) %>% 
  filter(!(is.na(count)) & count >= 5 & lat >= 25.5) %>% # & count > 50 
  left_join(weather.df %>% mutate(date = as.Date(date)) %>% 
              select(date,lat,lon,year,tailwind,precip7,temp7,windspeed7,winddir7),by = c("date","lat","lon","year")) %>%  
  left_join(ndvi.df %>% 
              select(julian,lat,lon,year,ndvi), by = c("julian","lat","lon","year")) %>% 
  mutate(std_yr = year - min(year),
         obs = row_number(),
         year_idx = as.numeric(factor(year)),
         site_idx = as.numeric(factor(paste(lat,lon))),
         temp = scale_this(temp7),
         rain = scale_this(precip7),
         tailwind_raw = tailwind,
         tailwind = scale_this(tailwind),
         rain_temp = scale_this(precip7*temp7),
         ndvi_raw = ndvi,
         ndvi = scale_this(ndvi),
         windspeed = scale_this(windspeed7),
         winddir = scale_this(winddir7),
         speed_dir = scale_this(windspeed7*winddir7))


# data filtering summary --------------------------------------------------
count_dat %>% 
  group_by(year) %>% 
  reframe(n = n()) %>% 
  pull(n) %>% 
  quantile(probs = c(0,.5,1))

# prepare model mesh ------------------------------------------------------
# make a set of distinct study sites for mapping
site_map <- count_dat %>%
  select(easting, northing) %>%
  distinct()

# make a two extension hulls and mesh for spatial model
hull <- fm_extensions(
  count_dat,
  convex = c(200, 500),
  concave = c(350, 500)
)
mesh <- fm_mesh_2d_inla(
  boundary = hull, max.edge = c(100, 600), # km inside and outside
  cutoff = 50, offset = c(100, 300),
  crs = fm_crs(count_dat)
) # cutoff is min edge
# save(mesh,file = "mesh.RData")
# plot the mesh
ggplot() +
  gg(data = mesh) +
  geom_sf(data = site_map, col = "darkgreen", size = 1) +
  geom_sf(data = map2 %>% st_simplify(dTolerance = 10,preserveTopology = T), fill = NA) +
  theme_bw() +
  labs(x = "", y = "")
# ggsave("plots/model mesh.png", bg = "white")


# prepare spde - spatial surface ------------------------------------------
# make spde
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(500, 0.5),
  prior.sigma = c(1, 0.5)
)
# save(spde,file = "spde.RData")

# additional priors ------------------------------------------------------------------
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))
# year svc model ----------------------------------------------------------
svc_formula <- count ~ .
svc_components <- ~ -1 +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+  
  tau(geometry, weights = std_yr, model = spde)

res.full <- bru(
  svc_components,
  like(
    formula = svc_formula,
    family = "nbinomial",
    data = count_dat
  ),
  options = list(
    control.predictor = list(compute = TRUE,quantiles = c(0.025,0.25,0.75,0.975)),
    control.compute = list(waic = TRUE, cpo = FALSE),
    control.inla = list(int.strategy = "eb"),
    verbose = FALSE
  )
)

res.full
res.full$summary.fixed
res.full$waic$waic
# save(res.full,file = "space time model results.RData")

# plot trend map ----------------------------------------------------------
# get easting and northing limits
bbox <- fm_bbox(hull[[1]])
grd_dims <- round(c(x = diff(bbox[[1]]), y = diff(bbox[[2]])) / 10)

# make mesh projector to get model summaries from the mesh to the mapping grid
mesh_proj <- fm_evaluator(
  mesh,
  xlim = bbox[[1]], ylim = bbox[[2]], dims = grd_dims
)

# pull data
kappa <- data.frame(
  median = exp(mod.res[[3]]$summary.random$kappa$"0.5quant"),
  range95 = exp(mod.res[[3]]$summary.random$kappa$"0.975quant") -
    exp(mod.res[[3]]$summary.random$kappa$"0.025quant")
)
alph <- data.frame(
  median = exp(mod.res[[3]]$summary.random$alpha$"0.5quant"),
  range95 = exp(mod.res[[3]]$summary.random$alpha$"0.975quant") -
    exp(mod.res[[3]]$summary.random$alpha$"0.025quant")
)
taus <- data.frame(
  median = (exp(mod.res[[3]]$summary.random$tau$"0.5quant") - 1) * 100,
  range95 = (exp(mod.res[[3]]$summary.random$tau$"0.975quant") -
               exp(mod.res[[3]]$summary.random$tau$"0.025quant")) * 100,
  ucl = (exp(mod.res[[3]]$summary.random$tau$"0.975quant") - 1) * 100,
  lcl = (exp(mod.res[[3]]$summary.random$tau$"0.025quant") - 1) * 100) %>% 
  mutate(med_dif0 = ifelse(ucl < 0 | lcl > 0,median,NA))

# loop to get estimates on a mapping grid
pred_grids <- lapply(
  list(alpha = alph, tau = taus),
  function(x) as.matrix(fm_evaluate(mesh_proj, x))
)


# make a terra raster stack with the posterior median and range95
out_stk0 <- rast()
for (j in 1:1) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                  Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j, crs = epsg6703km)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
                     Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j, crs = epsg6703km)
  out_j <- c(mean_j, range95_j)
  terra::add(out_stk0) <- out_j
}

for (j in 2:2) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                  Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j, crs = epsg6703km)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
                     Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j, crs = epsg6703km)
  mean_dif0_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                       Z = c(matrix(pred_grids[[j]][, 5], grd_dims[1]))
  )
  mean_dif0_j <- rast(mean_dif0_j, crs = epsg6703km)
  out_j <- c(mean_j, range95_j,mean_dif0_j)
  terra::add(out_stk0) <- out_j
}

names(out_stk0) <- c(
  "alpha_median", "alpha_range95",
  "tau_median", "tau_range95","tau_dif0"
)
out_stk0 <- terra::mask(out_stk0, map3, touches = FALSE)

make_plot_field <- function(data_stk, scale_label) {
  ggplot(map2 %>% st_simplify(dTolerance = 10,preserveTopology = T)) +
    geom_sf(fill = NA) +
    coord_sf(datum = NA) +
    geom_spatraster(data = data_stk) +
    labs(x = "", y = "") +
    scale_fill_distiller(scale_label,
                         palette = "Spectral",
                         na.value = "transparent",
                         direction = -1
    ) +
    # scale_fill_viridis_c(scale_label,na.value = "transparent")+
    geom_sf(data = site_map,size = .5,shape = 4,alpha = .5)+
    theme_bw()
}

make_plot_site <- function(data, scale_label) {
  ggplot(map2 %>% st_simplify(dTolerance = 10,preserveTopology = T)) +
    geom_sf() +
    coord_sf(datum = NA) +
    geom_sf(data = data, size = 1, mapping = aes(colour = value)) +
    scale_colour_distiller(scale_label, palette = "Spectral",direction = 1) +
    # scale_fill_viridis_c(scale_label,na.value = "transparent")+
    labs(x = "", y = "") +
    theme_bw() +
    geom_sf(fill = NA)
}
# medians
# fields alpha_s, tau_s
pa <- make_plot_field(
  data_stk = out_stk0[["alpha_median"]],
  scale_label = "posterior\nmedian\nexp(alpha_s)"
)
pt <- make_plot_field(
  data_stk = out_stk0[["tau_median"]],
  scale_label = "posterior\nmedian\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$median)),
  scale_label = "posterior\nmedian\nexp(kappa_s)"
)

# range95
# fields alpha_s, tau_s
pa_range95 <- make_plot_field(
  data_stk = out_stk0[["alpha_range95"]],
  scale_label = "posterior\nrange95\nexp(alpha_s)"
)
pt_range95 <- make_plot_field(
  data_stk = out_stk0[["tau_range95"]],
  scale_label = "posterior\nrange95\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps_range95 <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$range95)),
  scale_label = "posterior\nrange95\nexp(kappa_s)"
)

multiplot(ps,pa,pt, cols = 3)
multiplot(ps_range95,pa_range95,pt_range95, cols = 3)

lat.df0 = data.frame(taus = ((exp(mod.res2[[3]]$summary.random$tau$'0.5quant')-1)*100),
                     taus_ucl = ((exp(mod.res2[[3]]$summary.random$tau$'0.975quant')-1)*100),
                     taus_lcl = ((exp(mod.res2[[3]]$summary.random$tau$'0.025quant')-1)*100),
                     easting = mesh$loc[,1],
                     northing = mesh$loc[,2]) %>% 
  st_as_sf(coords = c("easting", "northing"), crs = epsg6703km, remove = FALSE) %>% 
  st_crop(st_bbox(out_stk0))

within10 <- sapply(st_within(lat.df0, map3), function(z) if (length(z)==0) 
  NA_integer_ else z[1])
lat.df02 = lat.df0 %>% 
  mutate(state = map3$ST_12[within10]) %>% 
  filter(!(is.na(state))) %>% 
  st_transform(crs = 4326) %>% 
  mutate(latitude = st_coordinates(.)[,2],
         longitude = st_coordinates(.)[,1])

p.lat0 = ggplot(lat.df02, aes(y = taus/100, x = latitude, color = longitude))+
  geom_point()+
  geom_smooth(method = "lm",se = F,color = "black")+
  ggpubr::theme_pubclean()+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  labs(y = "Estimated Annual Trend", x = "Latitude")+
  scale_color_viridis_c("Longitude")
p.lat0

hull.sig0 = lat.df0 %>% 
  mutate(state = map3$ST_12[within10]) %>% 
  filter(!(is.na(state))) %>% 
  mutate(taus_dif0 = ifelse(taus_ucl < 0 | taus_lcl > 0,taus,NA)) %>% 
  filter(taus > 0 | taus < 0) %>% 
  na.omit() %>% 
  st_geometry() %>% 
  st_combine() %>% 
  st_concave_hull(.05) %>% 
  st_as_sf()

p.tau0 = ggplot(map2 %>% st_simplify(dTolerance = 10,preserveTopology = T)) +
  geom_sf(fill = NA) +
  coord_sf(datum = NA) +
  geom_spatraster(data = out_stk0[["tau_median"]]/100) +
  # geom_sf(data = point.neg,fill = "black",alpha = 0.75, shape = 18) +
  geom_sf(data = hull.sig0,fill = NA,color = "black",size= 1.5) +
  labs(x = "", y = "") +
  # scale_fill_viridis_c(expression(paste(Delta," per year")),
  #                      na.value = "transparent",
  #                      labels = scales::label_percent())+
  scale_fill_distiller(expression(paste(Delta," per year")),
                       palette = "Spectral",
                       na.value = "transparent",
                       direction = 1,
                       labels = scales::label_percent()) +
  annotation_scale(location = "br", width_hint = 0.3) +
  geom_sf(inherit.aes = F, data = site_map,alpha=0.10,size = 1.25,pch=4)+
  theme_bw()+
  theme(legend.position = c(0.9, 0.45), legend.direction = "vertical",legend.background =element_rect(fill = "transparent"))
p.tau0

ggarrange(p.tau0,p.lat0, labels = c("(A)","(B)"))
ggsave("plots/monarch roost trend map and latitude.png",bg = "white",height = 4,width = 10)

p.tau0
ggsave("plots/monarch roost trend map and latitude.png",bg = "white",height = 4,width = 10)

# model selection -------------------------------------------------------
source("fixed model list.R")
model_list
model_names
mod.res = list()
for(i in 1:length(model_names)){
  svc_formula <- count ~ .
  mod.res[[i]] <- bru(
    model_list[[i]],
    like(
      formula = svc_formula,
      family = "nbinomial",
      data = count_dat
    ),
    options = list(
      control.predictor = list(compute = TRUE,quantiles = c(0.025,0.25,0.5,0.75,0.975)),
      control.compute = list(waic = TRUE, cpo = FALSE),
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE
    )
  )
}
waic.tab <- data.frame(mod_num = 1:length(model_names),model = model_names,waic = 1:length(model_names))
for(i in 1:nrow(waic.tab)){
  waic.tab$waic[i]<-ifelse(is.null(mod.res[[i]]$waic),NA,mod.res[[i]]$waic$waic[[1]])
}

waic.tab$dwaic = waic.tab$waic-min(waic.tab$waic)
waic.tab$rel_lkhood = exp(waic.tab$dwaic*-0.5)
waic.tab$weight = waic.tab$rel_lkhood/sum(waic.tab$rel_lkhood)
write_csv(waic.tab %>% mutate(waic = round(waic,2),
                              dwaic = round(dwaic,2),
                              rel_lkhood = round(rel_lkhood,3),
                              weight = round(weight,3)) %>% 
            arrange(dwaic),"final model selection results_space time drivers.csv")
save(mod.res,file = "output_all models_model selection.RData")

# model average -----------------------------------------------------------
#summarize model selection results
# load("output_all models_model selection.RData")
var_comp = c("temp","rain","wind","ndvi","rain*temp")
dwaic_thresh = waic.tab %>% filter(mod_num == 3) %>% pull(dwaic)
waic.tab.sum = data.frame(variable = var_comp,perc_present = 1:length(var_comp))
for(i in 1:length(var_comp)){
  mod_present = waic.tab %>% 
    filter(dwaic < dwaic_thresh) %>% 
    pull(model) %>%
    grep(.,pattern = ifelse(i == 5,"rain\\*temp",var_comp[i])) %>% 
    length()
  mod_total = waic.tab %>% 
    filter(dwaic < dwaic_thresh) %>% 
    pull(model) %>% 
    length()
  waic.tab.sum$perc_present[i] = (mod_present/mod_total)*100
}

#first, obtain list of models that contain a given variable and then get the position of that variable within the model that it is present in (index). I need the index because the inla output for the posterior marginals (as implemented via the bru()) gives posterior marginals labeled as index numbers that correspond to the position of the fixed effect variable in the formula. 
#update variables to examine
var_comp2 = waic.tab.sum$variable[which(waic.tab.sum$perc_present > 0)]
var_comp2[which(var_comp2 == "wind")] = "tailwind"
waic.tab = waic.tab %>% 
  mutate(model = str_replace(waic.tab$model,"wind","tailwind"))
mods.keep = waic.tab %>% filter(dwaic < dwaic_thresh) %>% pull(mod_num)
waic.tab_filt = waic.tab %>% 
  filter(dwaic < dwaic_thresh) # only average models that outperform intercept-only model
mods.to.average = list()
for(i in 1:length(var_comp2)){
  mods.to.average[[i]] = data.frame(var_present = waic.tab_filt$mod_num[waic.tab_filt$model %>% 
                                                                          grep(.,pattern = ifelse(i == 5,"rain\\*temp",var_comp[i]))])
  var_index = vector()
  for(j in 1:length(mods.to.average[[i]]$var_present)){
    mod.filt = mods.to.average[[i]]$var_present[j]
    var_index[j] = mod.res[[mod.filt]]$summary.fixed %>% 
      rownames() %>% 
      grep(.,pattern = ifelse(i == 5,"rain_temp" ,paste0("^",var_comp2[i],"$")))
  }
  mods.to.average[[i]]$index = var_index
  names(mods.to.average)[i] = var_comp2[i]
}

post.dist.new = list()
for(i in 1:length(var_comp2)){
  posterior.dens = list()
  for(j in 1:nrow(waic.tab_filt)){
    mod_num = waic.tab_filt$mod_num[j]
    mod = mod.res[[mod_num]]
    if(var_comp2[i] %in% (mod$summary.fixed %>% rownames())){
      index_num = mods.to.average[[var_comp2[i]]]$index[which(mods.to.average[[var_comp2[i]]]$var_present == mod_num)]
      mod_weight = (waic.tab_filt$weight[which(waic.tab_filt$mod_num == mod_num)])/sum(waic.tab_filt$weight)
      # x_vals = seq(from = min(xx), to = max(xx),length.out = 100)
      x_vals = seq(from = -3,to = 3,length.out = 50000)
      y_vals = mod_weight*(inla.dmarginal(x = x_vals,marginal = mod$marginals.fixed[[index_num]]))
      posterior.dens[[j]] = data.frame(mod_num = mod_num,
                                       mod_weight = mod_weight,
                                       x = x_vals,
                                       y = y_vals)
    } else{
      mod_weight = (waic.tab_filt$weight[which(waic.tab_filt$mod_num == mod_num)])/sum(waic.tab_filt$weight)
      # x_vals = seq(from = min(xx), to = max(xx),length.out = 100)
      x_vals = seq(from = -3,to = 3,length.out = 50000)
      y_vals = dnorm(x = x_vals, mean = 0, sd = 1)*mod_weight
      posterior.dens[[j]] = data.frame(mod_num = mod_num,
                                       mod_weight = mod_weight,
                                       x = x_vals,
                                       y = y_vals)
    }
  }
  post.dist.new[[var_comp2[i]]] = do.call("rbind",posterior.dens) %>% 
    group_by(x) %>% 
    reframe(y_new = sum(y)) %>% 
    select(y_new,x)
}

#model averaged covariate effects
##get vector beta estimates from marginals
samp.list = list()
for(i in 1:length(var_comp2)){
  dist.df = post.dist.new[[var_comp2[i]]] %>% 
    mutate(y = y_new*10)
  yy = c()
  for(j in 1:50000){
    yy = c(yy,rep(dist.df$x[j],dist.df$y[j]))
  }
  samp.list[[var_comp2[i]]] = sample(yy,5000,replace = TRUE) 
}

#histogram of averages posterior marginal distributions for fixed effect
hist.list = list()
for(i in 1:length(var_comp2)){
  sampx = samp.list[[var_comp2[i]]]
  sampmin = quantile(samp.list[[var_comp2[i]]],prob = c(0.025))
  sampmax = quantile(samp.list[[var_comp2[i]]],prob = c(0.975))
  perc_dif = ifelse(mean(sampx) > 0, 
                    (sum(sampx > 0)/length(sampx))*100, 
                    (sum(sampx < 0)/length(sampx))*100) %>% 
    round(., digits = 1)
  #I am trimming the samples for visualization purposes--the long tails due to shrinkage towards 0 makes it hard to see the marginal posteriors
  hist.list[[var_comp2[i]]] = ggplot(post.dist.new[[var_comp2[i]]] %>% filter(x > sampmin & x < sampmax),aes(y = y_new,x = x))+
    geom_line(linewidth = 1.5)+
    theme_pubclean()+
    geom_vline(xintercept = 0, color = "red", linewidth = 1.5)+
    labs(y = bquote("Pr(" ~ beta ~ ")"), x = bquote(beta), title = paste0("variable: ", var_comp2[i]))+
    annotate(geom = "text", x = Inf, y = Inf, hjust = "inward", vjust = "inward", label = paste0(perc_dif, ifelse(mean(samp.list[[var_comp2[i]]]) > 0,"% > 0", "% < 0")))
}
names(hist.list)
ggarrange(hist.list[["temp"]],hist.list[["rain"]],hist.list[["rain*temp"]],hist.list[["tailwind"]],hist.list[["ndvi"]], ncol = 3,nrow = 2, align = "hv",labels = c("A","B","C","D","E"))
ggsave("plots/covariate effects_weighted averages of marginal posteriors.png", 
       bg = "white",
       height = 7,
       width = 10)

# plot covariate effects using weighted averages for estimates ------------
#temp
p.temp = data.frame(temp = seq(from = min(count_dat$temp,na.rm = T),to = max(count_dat$temp,na.rm = T),length.out = 100)) %>%
  tidyr::expand(temp)  %>% 
  rownames_to_column() %>%   
  slice(rep(1:n(), each = 5000)) %>%
  mutate(b.temp= rep(samp.list[["temp"]],100),
         count = (b.temp*temp)) %>%
  group_by(temp) %>%
  reframe(count_med = median(count),count_ucl = quantile(count,probs = 0.75),count_lcl = quantile(count,probs = 0.25)) %>% 
  mutate(temp = (((temp*sd(count_dat$temp7,na.rm = T))+mean(count_dat$temp7,na.rm = T))-32) * (5/9)) %>% 
  ggplot(aes(y = count_med,x = temp))+
  geom_ribbon(aes(ymin = count_lcl,ymax = count_ucl),alpha = .25,linewidth = .1)+
  theme_pubclean()+
  geom_line(lty = "solid",linewidth = 1.15)+
  labs(y = "Relative Roost Size",
       x = "Temperature (°C)")
p.temp

#ndvi
p.ndvi = data.frame(ndvi = seq(from = min(count_dat$ndvi,na.rm = T),to = max(count_dat$ndvi,na.rm = T),length.out = 100)) %>% 
  slice(rep(1:n(), each = 5000)) %>% 
  rownames_to_column()  %>% 
  mutate(b.ndvi= rep(samp.list[["ndvi"]],100),
         count = ndvi*b.ndvi) %>% 
  group_by(ndvi) %>%
  reframe(count_med = median(count),count_ucl = quantile(count,probs = 0.75),count_lcl = quantile(count,probs = 0.25)) %>% 
  mutate(ndvi = (ndvi*sd(count_dat$ndvi_raw,na.rm = T))+mean(count_dat$ndvi_raw,na.rm = T)) %>% 
  ggplot(aes(y = count_med,x = ndvi))+
  geom_ribbon(aes(ymin = count_lcl,ymax = count_ucl),alpha = .25,linewidth = .1)+
  theme_pubclean()+
  geom_line(lty = "solid",linewidth = 1.15)+
  labs(y = "", 
       x = "NDVI")
p.ndvi

ggarrange(p.temp,p.ndvi,ncol = 2)
ggsave("plots/covariate effects.png",bg = "white",width = 7,height = 3)

# changes in covariates over time -----------------------------------------
library(glmmTMB)
library(sjPlot)
library(car)
str(count_dat)
mod.ndvi = glmmTMB(ndvi_raw ~ lat*std_yr,data = count_dat)
Anova(mod.ndvi)
summary(mod.ndvi)
ndvi.trend = plot_model(mod.ndvi,type = "pred",terms = c("std_yr","lat"))+
  labs(title = NULL, x = "Year (0 = 2007)", y = "NDVI", color = "Latitude")+
  theme_pubclean()

mod.temp = glmmTMB(temp ~ lat*std_yr,data = count_dat %>% 
                     mutate(temp = (temp7-32) * (5/9)))
Anova(mod.temp)
summary(mod.temp)
temp.trend = plot_model(mod.temp,type = "pred",terms = c("std_yr","lat"))+
  labs(title = NULL, x = "Year (0 = 2007)", y = "Temperature (°C)", color = "Latitude")+
  theme_pubclean()

ggarrange(temp.trend, ndvi.trend, common.legend = T)
ggsave("plots/temperature and ndvi trends.png",bg = "white",height = 4,width = 10)

# figure for phenology change over time -----------------------------------
mod.pheno = glmmTMB(julian~lat+std_yr, count_dat)
summary(mod.pheno)
car::Anova(mod.pheno)
plot_model(mod.pheno,type = "pred",terms = c("std_yr","lat"))+
  theme_pubclean()+
  labs(title = NULL, y = "Julian Day", x = "Year (0 = 2007)", color = "Latitude")
ggsave("plots/julian day vs year.png",bg = "white",height = 4,width = 4.5)


# view top model ----------------------------------------------------------
# load("fixed effect model selection results list.RData")
mod.res2[[30]]$summary.hyperpar %>% 
  select(mean,sd) %>% 
  rownames_to_column("parameter") %>% 
  mutate(mean = round(mean,3),
         sd = round(sd,3)) %>% 
  write_csv("range parameters table.csv",)

mod.res2[[30]]$summary.hyperpar

median(beta.df %>% filter(params == "Temperature") %>% pull(betas))/abs(mod.res2[[30]]$summary.random$tau$'0.5quant' %>% median())
median(beta.df %>% filter(params == "NDVI") %>% pull(betas))/abs(mod.res2[[30]]$summary.random$tau$'0.5quant' %>% median())
median(beta.df %>% filter(params == "Temperature") %>% pull(betas))/median(beta.df %>% filter(params == "NDVI") %>% pull(betas))
((exp(mod.res2[[30]]$summary.random$tau$'0.5quant')-1)*100) %>% quantile(prob = c(0.025,0.5,0.975))
exp(mod.res2[[30]]$summary.random$tau$'0.5quant')%>% quantile(prob = c(0.025,0.5,0.975))
(1-(0.9340051^17))*100
(1-(0.8873992^17))*100

count_dat %>% 
  mutate(countx = mod.res2[[30]]$summary.linear.predictor$mean[1:nrow(count_dat)]) %>% 
  ggplot(aes(y= countx,x= year, col = lat)) +
  geom_point()+
  geom_smooth(method = "lm")
# plot trend map from top model ----------------------------------------------------------
# get easting and northing limits
bbox <- fm_bbox(hull[[1]])
grd_dims <- round(c(x = diff(bbox[[1]]), y = diff(bbox[[2]])) / 10)

# make mesh projector to get model summaries from the mesh to the mapping grid
mesh_proj <- fm_evaluator(
  mesh,
  xlim = bbox[[1]], ylim = bbox[[2]], dims = grd_dims
)

# pull data
kappa <- data.frame(
  median = exp(mod.res[[30]]$summary.random$kappa$"0.5quant"),
  range95 = exp(mod.res[[30]]$summary.random$kappa$"0.975quant") -
    exp(mod.res[[30]]$summary.random$kappa$"0.025quant")
)
alph <- data.frame(
  median = exp(mod.res[[30]]$summary.random$alpha$"0.5quant"),
  range95 = exp(mod.res[[30]]$summary.random$alpha$"0.975quant") -
    exp(mod.res[[30]]$summary.random$alpha$"0.025quant")
)
taus <- data.frame(
  median = (exp(mod.res[[30]]$summary.random$tau$"0.5quant") - 1) * 100,
  range95 = (exp(mod.res[[30]]$summary.random$tau$"0.975quant") -
               exp(mod.res[[30]]$summary.random$tau$"0.025quant")) * 100,
  ucl = (exp(mod.res[[30]]$summary.random$tau$"0.975quant") - 1) * 100,
  lcl = (exp(mod.res[[30]]$summary.random$tau$"0.025quant") - 1) * 100) %>% 
  mutate(med_dif0 = ifelse(ucl < 0 | lcl > 0,median,NA))

# loop to get estimates on a mapping grid
pred_grids <- lapply(
  list(alpha = alph, tau = taus),
  function(x) as.matrix(fm_evaluate(mesh_proj, x))
)


# make a terra raster stack with the posterior median and range95
out_stk <- rast()
for (j in 1:1) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                  Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j, crs = epsg6703km)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
                     Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j, crs = epsg6703km)
  out_j <- c(mean_j, range95_j)
  terra::add(out_stk) <- out_j
}

for (j in 2:2) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                  Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j, crs = epsg6703km)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
                     Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j, crs = epsg6703km)
  mean_dif0_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                       Z = c(matrix(pred_grids[[j]][, 5], grd_dims[1]))
  )
  mean_dif0_j <- rast(mean_dif0_j, crs = epsg6703km)
  out_j <- c(mean_j, range95_j,mean_dif0_j)
  terra::add(out_stk) <- out_j
}

names(out_stk) <- c(
  "alpha_median", "alpha_range95",
  "tau_median", "tau_range95","tau_dif0"
)
out_stk <- terra::mask(out_stk, map3, touches = FALSE)

make_plot_field <- function(data_stk, scale_label) {
  ggplot(map2 %>% st_simplify(dTolerance = 10,preserveTopology = T)) +
    geom_sf(fill = NA) +
    coord_sf(datum = NA) +
    geom_spatraster(data = data_stk) +
    labs(x = "", y = "") +
    scale_fill_distiller(scale_label,
                         palette = "Spectral",
                         na.value = "transparent",
                         direction = -1
    ) +
    # scale_fill_viridis_c(scale_label,na.value = "transparent")+
    geom_sf(data = site_map,size = .5,shape = 4,alpha = .5)+
    theme_bw()
}

make_plot_site <- function(data, scale_label) {
  ggplot(map2 %>% st_simplify(dTolerance = 10,preserveTopology = T)) +
    geom_sf() +
    coord_sf(datum = NA) +
    geom_sf(data = data, size = 1, mapping = aes(colour = value)) +
    scale_colour_distiller(scale_label, palette = "Spectral",direction = 1) +
    # scale_fill_viridis_c(scale_label,na.value = "transparent")+
    labs(x = "", y = "") +
    theme_bw() +
    geom_sf(fill = NA)
}
# medians
# fields alpha_s, tau_s
pa <- make_plot_field(
  data_stk = out_stk[["alpha_median"]],
  scale_label = "posterior\nmedian\nexp(alpha_s)"
)
pt <- make_plot_field(
  data_stk = out_stk[["tau_median"]],
  scale_label = "posterior\nmedian\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$median)),
  scale_label = "posterior\nmedian\nexp(kappa_s)"
)

# range95
# fields alpha_s, tau_s
pa_range95 <- make_plot_field(
  data_stk = out_stk[["alpha_range95"]],
  scale_label = "posterior\nrange95\nexp(alpha_s)"
)
pt_range95 <- make_plot_field(
  data_stk = out_stk[["tau_range95"]],
  scale_label = "posterior\nrange95\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps_range95 <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$range95)),
  scale_label = "posterior\nrange95\nexp(kappa_s)"
)

png("plots/spatial fields of paramaters.png",pointsize =16, height = 3000, width = 4500,res = 300)
multiplot(ps,ps_range95,pa,pa_range95,pt,pt_range95, cols = 3)
dev.off()
#figure 1
lat.df = data.frame(taus = ((exp(mod.res[[30]]$summary.random$tau$'0.5quant')-1)*100),
                    taus_ucl = ((exp(mod.res[[30]]$summary.random$tau$'0.975quant')-1)*100),
                    taus_lcl = ((exp(mod.res[[30]]$summary.random$tau$'0.025quant')-1)*100),
                    easting = mesh$loc[,1],
                    northing = mesh$loc[,2]) %>% 
  st_as_sf(coords = c("easting", "northing"), crs = epsg6703km, remove = FALSE) %>% 
  st_crop(st_bbox(out_stk))

within1 <- sapply(st_within(lat.df, map3), function(z) if (length(z)==0) 
  NA_integer_ else z[1])
lat.df2 = lat.df %>% 
  mutate(state = map3$ST_12[within1]) %>% 
  filter(!(is.na(state))) %>% 
  st_transform(crs = 4326) %>% 
  mutate(latitude = st_coordinates(.)[,2],
         longitude = st_coordinates(.)[,1])

p.lat = ggplot(lat.df2, aes(y = taus/100, x = latitude, color = longitude))+
  geom_point()+
  geom_smooth(method = "lm",se = F,color = "black")+
  ggpubr::theme_pubclean()+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  labs(y = "Estimated Annual Trend", x = "Latitude")+
  scale_color_viridis_c("Longitude")
p.lat

hull.sig = lat.df %>% 
  mutate(state = map3$ST_12[within1]) %>% 
  filter(!(is.na(state))) %>% 
  mutate(taus_dif0 = ifelse(taus_ucl < 0 | taus_lcl > 0,taus,NA)) %>% 
  filter(taus > 0 | taus < 0) %>% 
  na.omit() %>% 
  st_geometry() %>% 
  st_combine() %>% 
  st_concave_hull(.05) %>% 
  st_as_sf()

p.tau = ggplot(map2 %>% st_simplify(dTolerance = 10,preserveTopology = T)) +
  geom_sf(fill = NA) +
  coord_sf(datum = NA) +
  geom_spatraster(data = out_stk[["tau_median"]]/100) +
  # geom_sf(data = point.neg,fill = "black",alpha = 0.75, shape = 18) +
  geom_sf(data = hull.sig,fill = NA,color = "black",size= 1.5) +
  labs(x = "", y = "") +
  # scale_fill_viridis_c(expression(paste(Delta," per year")),
  #                      na.value = "transparent",
  #                      labels = scales::label_percent())+
  scale_fill_distiller(expression(paste(Delta," per year")),
                       palette = "Spectral",
                       na.value = "transparent",
                       direction = 1,
                       labels = scales::label_percent()) +
  annotation_scale(location = "br", width_hint = 0.3) +
  geom_sf(inherit.aes = F, data = site_map,alpha=0.10,size = 1.25,pch=4)+
  theme_bw()+
  theme(legend.position = c(0.9, 0.45), legend.direction = "vertical",legend.background =element_rect(fill = "transparent"))
p.tau

ggarrange(p.tau,p.lat, labels = c("(A)","(B)"))
ggsave("plots/monarch roost trend map and latitude_with covariates.png",bg = "white",height = 4,width = 10)

p.tau
ggsave("plots/monarch roost trend map_with covariates.png",bg = "white",height = 4.25,width = 5.25)


# compare tau with and without covariates ---------------------------------
ggarrange(p.tau0,p.tau)
ggsave("plots/monarch roost trend map_ without and with covariates.png",bg = "white",height = 4,width = 10)



