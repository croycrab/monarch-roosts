#model selection possible models
model_names = c("intercept","fixed_year","svc_year",
                "fixed_wind","fixed_rain","fixed_temp","fixed_ndvi",
                "fixed_wind_rain","fixed_wind_temp","fixed_wind_ndvi","fixed_rain_temp","fixed_rain_ndvi","fixed_temp_ndvi",
                "fixed_wind_rain_temp","fixed_wind_rain_ndvi","fixed_rain_temp_ndvi",
                "fixed_wind_rain_temp_ndvi",
                "fixed_wind_svc_year","fixed_rain_svc_year","fixed_temp_svc_year","fixed_ndvi_svc_year",
                "fixed_wind_rain_svc_year","fixed_wind_temp_svc_year","fixed_wind_ndvi_svc_year","fixed_rain_temp_svc_year","fixed_rain_ndvi_svc_year","fixed_temp_ndvi_svc_year",
                "fixed_wind_rain_temp_svc_year","fixed_wind_rain_ndvi_svc_year","fixed_rain_temp_ndvi_svc_year",
                "fixed_wind_rain_temp_ndvi_svc_year",
                "fixed_rain*temp", "fixed_rain*temp_wind","fixed_rain*temp_ndvi","fixed_rain*temp_wind_ndvi","fixed_rain*temp_svc_year", "fixed_rain*temp_wind_svc_year","fixed_rain*temp_ndvi_svc_year","fixed_rain*temp_wind_ndvi_svc_year"
)

model_list = list() 
model_list[[1]] = ~ -1 +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[2]] = ~ -1 + std_yr + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[3]] = ~ -1 +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)


# single fixed --------------------------------------------------------------
model_list[[4]] = ~ -1 + tailwind +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[5]] = ~ -1 + rain +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[6]] = ~ -1 + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[7]] = ~ -1 + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

# pairwise fixed effects -----------------------------------------------------------
model_list[[8]] = ~ -1 + tailwind + rain +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[9]] = ~ -1 + tailwind + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[10]] = ~ -1 + tailwind + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[11]] = ~ -1 + rain + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[12]] = ~ -1 + rain + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[13]] = ~ -1 + temp + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)


# three-way fixed effects ----------------------------------------------------------
model_list[[14]] = ~ -1 + tailwind + rain + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[15]] = ~ -1 + tailwind + rain + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[16]] = ~ -1 + rain + temp + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

# four-way fixeds -----------------------------------------------------------
model_list[[17]] = ~ -1 + tailwind + rain + temp + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

# single fixed effects + year ---------------------------------------------
model_list[[18]] = ~ -1 + tailwind +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[19]] = ~ -1 + rain +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[20]] = ~ -1 + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[21]] = ~ -1 + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

# pairwise fixed effects + year -----------------------------------------------------------
model_list[[22]] = ~ -1 + tailwind + rain +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[23]] = ~ -1 + tailwind + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[24]] = ~ -1 + tailwind + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[25]] = ~ -1 + rain + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[26]] = ~ -1 + rain + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[27]] = ~ -1 + temp + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

# three-way fixed effects + year ----------------------------------------------------------
model_list[[28]] = ~ -1 + tailwind + rain + temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[29]] = ~ -1 + tailwind + rain + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[30]] = ~ -1 + rain + temp + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

# four-way fixed effects + year -----------------------------------------------------------
model_list[[31]] = ~ -1 + tailwind + rain + temp + ndvi + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)



# rain by temp interaction ------------------------------------------------
model_list[[32]] = ~ -1 + rain + temp +  rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[33]] = ~ -1 + tailwind + rain + temp +  rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[34]] = ~ -1 + rain + temp + ndvi +  rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[35]] = ~ -1 + tailwind + rain + temp + ndvi +  rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)

model_list[[36]] = ~ -1 + rain + temp + rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[37]] = ~ -1 + tailwind + rain + temp +  rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[38]] = ~ -1 + rain + temp + ndvi +  rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)

model_list[[39]] = ~ -1 + tailwind + rain + temp + ndvi +  rain_temp + 
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde)+
  tau(geometry, weights = std_yr, model = spde)





