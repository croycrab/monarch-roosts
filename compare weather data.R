#compare weather data from different sources

#libraries
library("tidyverse")

#import data
prism.df = read_csv("prism.csv") %>% 
  rename(id = ID) %>%
  group_by(id) %>% 
  reframe(tmean_prism = mean(tmean, na.rm = T),
          ppt_prism = sum(ppt, na.rm = T)) %>% 
  na.omit()

daymet.df = read_csv("daymet.csv") %>% 
  mutate(tmean = (tmin + tmax)/2) %>% 
  group_by(id) %>% 
  reframe(tmean_daymet = mean(tmean, na.rm = T),
          ppt_daymet = sum(precip, na.rm = T)) %>% 
  na.omit()

vc.df = read_csv("curated weather.csv") %>% 
  rename(ppt_vc = precip7,
         tmean_vc = temp7) %>% 
  select(id, tmean_vc, ppt_vc)

weather.df = vc.df %>% 
  left_join(prism.df, by = c("id")) %>% 
  left_join(daymet.df, by = c("id"))

p.t1 = ggplot(weather.df, aes(x = (tmean_vc-32)*(5/9), y = tmean_daymet))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic() +
  labs(y = "Temperature (°C) -- DAYMET", x = "Temperature (°C) -- VC") +
  annotate("text",x = 10, y = Inf, label = paste0("r = ", cor(weather.df$tmean_vc,weather.df$tmean_daymet,use = 'complete.obs') %>% 
             round(.,2)), hjust = 0, vjust = 1)
# cor(weather.df$tmean_vc,weather.df$tmean_daymet,use = 'complete.obs')

p.t2 = ggplot(weather.df, aes(x = (tmean_vc-32)*(5/9), y = tmean_prism))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  labs(y = "Temperature (°C) -- PRISM", x = "Temperature (°C) -- VC") +
  annotate("text",x = 10, y = Inf, label = paste0("r = ", cor(weather.df$tmean_vc,weather.df$tmean_prism,use = 'complete.obs') %>% 
                                                    round(.,2)), hjust = 0, vjust = 1)
# cor(weather.df$tmean_vc,weather.df$tmean_prism,use = 'complete.obs')

p.p1 = ggplot(weather.df, aes(x = ppt_vc*25.2, y = ppt_daymet))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  labs(y = "Precipitation (mm) -- DAYMET", x = "Precipitation (mm) -- VC") +
  annotate("text",x = 25.4, y = Inf, label = paste0("r = ", cor(weather.df$ppt_vc,weather.df$ppt_daymet,use = 'complete.obs') %>% 
                                                    round(.,2)), hjust = 0, vjust = 1)
# cor(weather.df$ppt_vc,weather.df$ppt_daymet,use = 'complete.obs')

p.p2 = ggplot(weather.df, aes(x = ppt_vc*25.4, y = ppt_prism))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  labs(y = "Precipitation (mm) -- PRISM", x = "Precipitation (mm) -- VC") +
  annotate("text",x = 25.4, y = Inf, label = paste0("r = ", cor(weather.df$ppt_vc,weather.df$ppt_prism,use = 'complete.obs') %>% 
                                                    round(.,2)), hjust = 0, vjust = 1)
  # cor(weather.df$ppt_vc,weather.df$ppt_prism,use = 'complete.obs')

ggpubr::ggarrange(p.t1,p.t2,p.p1,p.p2,nrow = 2, ncol = 2, align = "hv")
ggsave("plots/weather comaprison.png", height = 6, width = 7)
