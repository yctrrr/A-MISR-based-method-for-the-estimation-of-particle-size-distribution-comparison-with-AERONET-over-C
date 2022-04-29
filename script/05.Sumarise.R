library(summarytools)
library(dplyr)

source("Spatial_analysis.R")

##  China polygon
polygon <- readOGR(dsn="data/other_data/prov_China/bou2_4p_region.shp",
                   layer="bou2_4p_region")
polygon@proj4string <-  CRS("+proj=longlat +datum=WGS84")
dir_predict <- "result/Table/MISR_AOD_Distribution_Prediction/"
##------01.Summary prediction results-------
daily <- fread(paste0(dir_predict,"Daily/MISR_particle_volume_",year,".csv"))

multiyear_mean_V_result <- fread(paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_mean.csv"))
multiyear_mean_V_ratio_result <- multiyear_mean_V_result%>%
  group_by(x,y) %>% 
  summarise(V_mean = V_mean[radius == 0.35]/V_mean[radius == 0.7])%>%
  mutate(radius = 0)

aod_multiyear_mean <- full_join(multiyear_mean_V_result, multiyear_mean_V_ratio_result)
aod_multiyear_mean
SR_file_name <- c("jingjinji","changsanjiao","zhusanjiao")
region_name <- c("BTH", "YRD", "PRD")

for (i in 1:length(SR_file_name)) {
  SR <- readOGR(paste0("data/other_data/specific_region/",SR_file_name[i],"/",SR_file_name[i],".shp"))
  assign(paste0("SR_",i),SR)
  multiyear_SR <- points_mask(aod_multiyear_mean, polygon = SR)%>%
    mutate(region = region_name[i])
  multiyear_SR
  mean(multiyear_SR$v_value)
  if(i == 1){
    multiyear_SR_join <- multiyear_SR
  }else{
    multiyear_SR_join <- full_join(multiyear_SR_join, multiyear_SR)
  }
}

multiyear_summary <- aod_multiyear_mean %>% 
  group_by(radius) %>% 
  descr("V_mean")%>%
  tb()

multiyear_summary
fwrite(multiyear_summary, file = paste0("result/Table/Summary/multiyear_summary.csv"))

##  group by economic regions
multiyear_SR_summary <- multiyear_SR_join %>%
  group_by(radius, region) %>% 
  descr("V_mean")%>%
  tb()
fwrite(multiyear_SR_summary, file  = paste0("result/Table/Summary/multiyear_regions_summary.csv"))

##  join province names
V_sf <- st_as_sf(aod_multiyear_mean, coords = c("x","y"), crs = "+proj=longlat +datum=WGS84")
polygon_sf <- st_as_sf(polygon)[,"NAME"]
V_sf_join <- st_join(V_sf, polygon_sf, join = st_nearest_feature, left = T)

##  group by provinces
multiyear_prov_summary <- V_sf_join %>%
  group_by(radius, NAME) %>% 
  descr("V_mean")%>%
  tb()
multiyear_prov_summary
fwrite(multiyear_prov_summary, file  = paste0("result/Table/Summary/multiyear_prov_summary.csv"))

##------Summary seasonal ratios-------
multiseason_mean_V_result <- readRDS(paste0(dir_predict, "Seasonal/MISR_GRID_particle_volume_multiseason_join.Rds"))
breaks <- c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.001)
multiseason_mean_V_ratio_result <- multiseason_mean_V_result%>%
  group_by(x, y, season, radius)%>%
  summarise(across("V_mean", mean))%>%
  summarise(V_ratio = V_mean[radius == 0.35]/V_mean[radius == 0.7])%>%
  mutate(discrete = cut(V_ratio, breaks = breaks, right = T, 
                        labels =  c("< 0.4","0.4 - 0.5","0.5 - 0.6","0.6 - 0.7","0.7 - 0.8","0.8 - 0.9","> 0.9")))

seasonal_ratio_summary <- multiseason_mean_V_ratio_result %>% 
  group_by(season) %>% 
  descr("V_ratio")%>%
  tb()

fwrite(seasonal_ratio_summary, file = paste0("result/Table/Summary/seasonal_ratio_summary.csv"))

##  join province names
ratio_sf <- st_as_sf(multiseason_mean_V_ratio_result, coords = c("x","y"), crs = "+proj=longlat +datum=WGS84")
polygon_sf <- st_as_sf(polygon)[,"NAME"]
ratio_sf_join <- st_join(ratio_sf, polygon_sf, join = st_nearest_feature, left = T)

seasonal_prov_summary <- ratio_sf_join %>%
  group_by(season, NAME) %>% 
  descr("V_ratio")%>%
  tb()
seasonal_prov_summary
fwrite(seasonal_prov_summary, file  = paste0("result/Table/Summary/seasonal_ratio_prov_summary.csv"))
