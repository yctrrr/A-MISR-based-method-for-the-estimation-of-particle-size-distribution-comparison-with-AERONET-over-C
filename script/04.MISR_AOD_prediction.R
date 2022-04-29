library(sp)
library(sf)
library(raster)
library(dplyr)
library(rgdal)
library(ggplot2)
library(data.table)
library(SearchTrees)
library(purrr)
library(Cairo)
library(colorspace)
library(grid)

source("AOD_to_VSD.R")
source("Scatter_plot.R")
source("Spatial_analysis.R")
source("Grid_visualization.R")
##  load in MISR AOD densitity
load("result/Rdata/MISR_AOD_Parameter/MISR_AOD_SD.Rdata")

##  China polygon
polygon <- readOGR(dsn="data/other_data/prov_China/bou2_4p_region.shp",
                   layer="bou2_4p_region")
polygon@proj4string <-  CRS("+proj=longlat +datum=WGS84")

integral <- function(density,lower,upper,n){
  sel <- density$x <= upper&density$x >= lower
  integral <- sum(frollmean(density$x[sel],2)^n*diff(density$x[sel])*frollmean(density$y[sel],2),na.rm=T)
  return(integral)
}
##  normalize function
standard <- function(x){(x-min(x))/(max(x)-min(x))}

dir_predict <- "result/Table/MISR_AOD_Distribution_Prediction/"
##-----01.Calculate all MISR AOD distribution and volume---------
rm_valid <- c(0.35,0.7)
for (year in 2004:2016){
  cat("----Year:",year,"-----\n")
  ##  read in observed MISR AOD data
  aod <- fread(paste0("data/AOD/MISR AOD/China_MISR_AOT_component_AGP_",year,".csv"))
  ##  Calculate total AOD
  aod <- mutate(aod, Component_all = Component_01 + Component_02 + Component_03 + Component_06 +
                  Component_08 + Component_14 + Component_19 + Component_21)%>%
    filter(Component_all!=0)
  

  aod_calcu <- misr_aod_to_distri(aod,rm_valid = rm_valid)%>%
    mutate(V = standard(V), N = standard(N))%>%
    dplyr::select(Year,Month,Day,DOY,Hour,Minute,MISRLon,MISRLat,radius,V)
  fwrite(aod_calcu,paste0(dir_predict,"Daily/MISR_particle_volume_",year,".csv"))
  saveRDS(aod_calcu,paste0(dir_predict,"Daily/MISR_particle_volume_",year,".Rds"))
}

##-----02.Calculate monthly MISR AOD distribution and volume---------
rm_valid <- c(0.35,0.7)
extent <- c(70,140,10,59)
xlim<- extent[2] - extent[1]
ylim <- extent[4] - extent[3]
raster_nation <- raster(ncol = xlim/0.2, nrow = ylim/0.2)
extent(raster_nation) <- extent
raster_nation  <- setValues(raster_nation ,1:ncell(raster_nation))

names(raster_nation) <- "GRID_ID"
T.month <- formatC(c(1:12), width = 2, format = "d", flag = "0")
for (year in 2004:2016) {
  # aod_calcu <- fread(paste0(dir_predict,"Daily/MISR_particle_volume_",year,".csv"))
  aod_calcu <- readRDS(paste0(dir_predict,"Daily/MISR_particle_volume_",year,".Rds"))
  if(!dir.exists(paste0(dir_predict, "Monthly/",year))){
    dir.create(paste0(dir_predict, "Monthly/",year))
  }
  for (month in 1:12) {
    cat("----Year: ",year,"----Month: ",month,"-----\n")
    aod_monthly_calcu <- aod_calcu%>%
      filter(Month == month)%>%
      group_by(MISRLon,MISRLat,radius)%>%
      summarise(across(c("V"), list(mean=mean)))
    fwrite(aod_monthly_calcu, paste0(dir_predict,"Monthly/",year,"/MISR_Point_particle_volume_",year,T.month[month],"_monthly_mean.csv"))

    aod_monthly_griddf <- lapply(rm_valid, function(i){points_to_grids(aod_monthly_calcu%>%filter(radius==i), x = "MISRLon", y = "MISRLat", polygon = polygon, raster = raster_nation)})%>%
      bind_rows()
    
    fwrite(aod_monthly_griddf, paste0(dir_predict, "Monthly/",year,"/MISR_GRID_particle_volume_",year,T.month[month],"_monthly_mean.csv"))
    
    aod_monthly_griddf$year <- year
    aod_monthly_griddf$month <- month
    if(year == 2004 & month == 1){
      aod_multimonth_join<- aod_monthly_griddf
    }else{
      aod_multimonth_join <- full_join(aod_multimonth_join ,aod_monthly_griddf)
    }
  }
}
fwrite(aod_multimonth_join, paste0(dir_predict, "Monthly/MISR_GRID_particle_volume_multimonth_join.csv"))
saveRDS(aod_multimonth_join, paste0(dir_predict, "Monthly/MISR_GRID_particle_volume_multimonth_join.Rds"))

##-----03.Calculate annual and multiyear MISR AOD distribution and volume-------
rm_valid <- c(0.35,0.7)
options(digits=10)
for (year in 2004:2016) {
  aod_annual_calcu  <- readRDS(paste0(dir_predict,"Daily/MISR_particle_volume_",year,".Rds"))%>%
    # fread(paste0(dir_predict,"Daily/MISR_particle_volume_",year,".csv"))%>%
    select(MISRLon, MISRLat, radius,V)%>%
    group_by(MISRLon, MISRLat, radius)%>%
    summarise(across(c("V"), list(mean=mean)))

  fwrite(aod_annual_calcu, paste0(dir_predict,"Annual/MISR_Point_particle_volume_",year,"_annual_mean.csv"))

  aod_annual_griddf <- lapply(rm_valid, function(i){points_to_grids(aod_annual_calcu%>%filter(radius==i), x = "MISRLon", y = "MISRLat", polygon = polygon, raster = raster_nation)})%>%
    bind_rows()
  # aod_annual_griddf <- points_to_grids(aod_annual_calcu, x = "MISRLon", y = "MISRLat", polygon = polygon, raster = raster_nation)
  fwrite(aod_annual_griddf, paste0(dir_predict,"Annual/MISR_GRID_particle_volume_",year,"_annual_mean.csv"))
  
  aod_annual_griddf$year <- year
  if(year == 2004){
    aod_multiyear_join<- aod_annual_griddf
  }else{
    aod_multiyear_join <- full_join(aod_multiyear_join ,aod_annual_griddf)
  }
}
fwrite(aod_multiyear_join, paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_join.csv"))
saveRDS(aod_multiyear_join, paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_join.Rds"))

aod_multiyear_mean <- aod_multiyear_join%>%
  group_by(x,y,radius)%>%
  summarise(across(c("V_mean"), mean))
fwrite(aod_multiyear_mean, paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_mean.csv"))

##-----04.Calculate seasonal MISR AOD distribution and volume---------
rm_valid <- c(0.35,0.7,50)
season_name <- c("Spring","Summer","Autumn","Winter")
for (year in 2004:2016) {
  ##  categorize season id
  aod_calcu <- fread(paste0(dir_predict,"Daily/MISR_particle_volume_",year,".csv"))%>%
    mutate(season_id = ifelse(Month%in%c(3,4,5),1, 
                              ifelse(Month%in%c(6,7,8),2,
                                     ifelse(Month%in%c(9,10,11),3,4))))
  if(!dir.exists(paste0(dir_predict, "Seasonal/",year))){
    dir.create(paste0(dir_predict, "Seasonal/",year))
  }
  for (season in 1:4) {
    cat("----Year: ",year,"----Season: ",season,"-----\n")
    aod_seasonal_calcu <- aod_calcu%>%
      filter(season_id == season)%>%
      group_by(MISRLon,MISRLat,radius)%>%
      summarise(across(c("V"), list(mean=mean)))
    fwrite(aod_seasonal_calcu, paste0(dir_predict,"Seasonal/",year,"/MISR_Point_particle_volume_",year,"_",season_name[season],"_mean.csv"))
    
    aod_seasonal_griddf <- lapply(rm_valid, function(i){points_to_grids(aod_seasonal_calcu%>%filter(radius==i), x = "MISRLon", y = "MISRLat", polygon = polygon, raster = raster_nation)})%>%
      bind_rows()
    fwrite(aod_seasonal_griddf, paste0(dir_predict, "Seasonal/",year,"/MISR_GRID_particle_volume_",year,"_",season_name[season],"_mean.csv"))

    aod_seasonal_griddf$year <- year
    aod_seasonal_griddf$season <- season_name[season]

    if(year == 2004 & season == 1){
      aod_multiseason_join <- aod_seasonal_griddf
    }else{
      aod_multiseason_join <- full_join(aod_multiseason_join ,aod_seasonal_griddf)
    }
  }
}
aod_multiseason_join <- mutate(aod_multiseason_join, across(season, factor,levels = season_name))
fwrite(aod_multiseason_join, paste0(dir_predict, "Seasonal/MISR_GRID_particle_volume_multiseason_join.csv"))
saveRDS(aod_multiseason_join, paste0(dir_predict, "Seasonal/MISR_GRID_particle_volume_multiseason_join.Rds"))

##-------Plot histogram distribution----------
aod_multiyear_join <- fread(paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_join.csv"))
aod_multiyear_mean <- fread(paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_mean.csv"))
aod_multiyear_join
h_v1 <- ggplot(aod_multiyear_join%>%filter(radius == 0.5), aes(V_mean))+
  geom_histogram()
ggsave(paste0("result/Fig/MISR_AOD_Volume_and_Distribution_Prediction/multiyear_volume_1um_histgoram.png"),h_v1)

h_v2.5 <- ggplot(aod_multiyear_join%>%filter(radius == 1.25), aes(V_mean))+
  geom_histogram()
ggsave(paste0("result/Fig/MISR_AOD_Volume_and_Distribution_Prediction/multiyear_volume_2.5um_histgoram.png"),h_v2.5)


multiyear_V_ratio_result <- aod_multiyear_join%>%
  group_by(x,y) %>% 
  summarise(V_ratio = V_mean[radius == 0.5]/V_mean[radius == 1.25])

multiyear_mean_V_ratio_result <- aod_multiyear_mean%>%
  group_by(x,y) %>% 
  summarise(V_ratio = V_mean[radius == 0.5]/V_mean[radius == 1.25])

h_v_ratio <- ggplot(multiyear_V_ratio_result, aes(V_ratio))+
  geom_histogram()
ggsave(paste0("result/Fig/MISR_AOD_Volume_Ratio/multiyear_volume_ratio_histgoram.png"),h_v_ratio)

h_v_mean_ratio <- ggplot(multiyear_mean_V_ratio_result, aes(V_ratio))+
  geom_histogram()
ggsave(paste0("result/Fig/MISR_AOD_Volume_Ratio/multiyear_mean_volume_ratio_histgoram.png"),h_v_mean_ratio)

##-----05.Annual volume distribution plot---------
dir_predict <- "result/Table/MISR_AOD_Distribution_Prediction/"
breaks <- c(0, 0.3, 0.4, 0.5, 0.6, 0.7, 1.001)
for (year in 2004:2016) {
  V_result <- fread(paste0(dir_predict,"Annual/MISR_GRID_particle_volume_",year,"_annual_mean.csv"))
  V_ratio_result <- V_result%>%
    group_by(x,y) %>%
    summarise(V_ratio = V_mean[radius == 0.35]/V_mean[radius == 0.7])

##--------05.1 Volume plot------
  colors <- RColorBrewer::brewer.pal(6,"YlGnBu")
  V_plot <- grid_facet_plot(data = V_result, label = year, value = "V_mean", facet = "radius")+
    facet_wrap("radius", nrow = 1, labeller = labeller(radius = c("0.35" = "Particle with radius less than 0.35 μm",
                                                               "0.7" = "Particle with radius less than 0.7 μm"),
                                                    groupwrap = label_wrap_gen(10)))+
    geom_polygon(data=polygon,aes(x = long, y = lat, group = group),
                 colour='black', fill='white',alpha=0.1)+
    scale_fill_gradientn(name=expression("Normalized volume"),
                       colours = colors,
                       guide = guide_colorbar(frame.colour = "black",frame.linewidth=2,
                                              title.position="right",title.hjust=0.5,title.vjust=1,ticks= F))

  V_plot_gt = ggplot_gtable(ggplot_build(V_plot))
  V_plot_gt$heights = 1.3*V_plot_gt$heights

  CairoPNG(paste0("result/Fig/04.MISR_AOD_Volume/volume_plot_",year,".png"),width = 2800,height = 1500,res=150)
  print(grid.draw(V_plot_gt))
  dev.off()
##--------05.2 Ratio plot------
  colors <- colorspace::sequential_hcl(6, palette = "agGrnYl")
  V_ratio_result <- mutate(V_ratio_result,
                           discrete = cut(V_ratio, breaks = breaks, 
                                          abels =  c("< 0.3","0.3 - 0.4","0.5 - 0.6","0.4 - 0.5","0.6 - 0.7","> 0.7")))

  ratio_plot <- grid_plot(data = V_ratio_result, label = year, value = "discrete")+
    geom_polygon(data=polygon,aes(x = long, y = lat, group = group),
                 colour='black', fill='white',alpha=0.1)+
    scale_fill_manual(name = expression("Volume ratio"),
                      values =  colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2, 
                                           title.position="right",title.hjust=0.5,title.vjust=1,
                                           byrow = TRUE, ticks= F))

  ratio_plot_gt = ggplot_gtable(ggplot_build(ratio_plot))
  ratio_plot_gt$heights = 1.2*ratio_plot_gt$heights
  
  CairoPNG(paste0("result/Fig/05.MISR_AOD_Volume_Ratio/volume_ratio_plot_",year,".png"),width = 2800,height = 1500,res=150)
  print(grid.draw(ratio_plot_gt))
  dev.off()
}

##-----05.3.Multiyear volume and ratio plot-------
multiyear_V_result <- fread(paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_join.csv"))

multiyear_Vsmall_result <- filter(multiyear_V_result, radius == 0.35)%>%
  mutate(V_mean = ifelse(V_mean > 0.5, 0.5, V_mean))

multiyear_Vmedium_result <- filter(multiyear_V_result, radius == 0.7)%>%
  mutate(V_mean = ifelse(V_mean > 0.5, 0.5, V_mean))

multiyear_V_ratio_result <- multiyear_V_result%>%
  group_by(x,y,year) %>% 
  summarise(V_ratio = V_mean[radius == 0.35]/V_mean[radius == 0.7])
multiyear_V_ratio_result
# multiyear_V_ratio_result%>%summary()
colors <- RColorBrewer::brewer.pal(6,"YlGnBu")
##  plot multiyear volume of particle with small radius
multiyear_Vsmall_result
multiyear_Vsmall_plot <- grid_facet_plot(multiyear_Vsmall_result, value = "V_mean", facet = "year",facet_text = 18)+
  facet_wrap("year", nrow = 3, labeller = labeller(groupwrap = label_wrap_gen(10)))+
  geom_polygon(data=polygon,aes(x = long, y = lat, group = group),
               colour='black', fill='white', alpha=0.1)+
  scale_fill_gradientn(name=expression("Normalized volume"),
                       colours = colors,
                       trans = "sqrt",
                       guide = guide_colorbar(frame.colour = "black",frame.linewidth=2,
                                              title.position="right", title.hjust=0.5, title.vjust=1,ticks= F))+
  theme(legend.title = element_text(colour = "black", size = rel(3.5)))+
  theme(strip.text = element_text(size = 35))

multiyear_Vsmall_plot <- add_coordinates(multiyear_Vsmall_plot, xmin = 70, xmax = 140, xstep = 15, ymin = 10, ymax = 60, ystep = 10, size = 1)

multiyear_Vsmall_plot_gt = ggplot_gtable(ggplot_build(multiyear_Vsmall_plot))
multiyear_Vsmall_plot_gt$heights = 1.3*multiyear_Vsmall_plot_gt$heights

CairoPNG(paste0("result/Fig/04.MISR_AOD_Volume/multiyear_volume_small_plot.png"),width = 1800,height = 700, res=60)
print(grid.draw(multiyear_Vsmall_plot_gt))
dev.off()

##  plot multiyear volume of particle with medium diameter
multiyear_Vmedium_plot <- grid_facet_plot(multiyear_Vmedium_result, value = "V_mean", facet = "year",facet_text = 18)+
  facet_wrap("year", nrow = 3, labeller = labeller(groupwrap = label_wrap_gen(20)))+
  geom_polygon(data=polygon,aes(x = long, y = lat, group = group),
               colour='black', fill='white',alpha=0.1)+
  scale_fill_gradientn(name=expression("Normalized volume"),
                       colours = colors,
                       trans = "sqrt",
                       guide = guide_colorbar(frame.colour = "black",frame.linewidth=2, 
                                              title.position="right",title.hjust=0.5,title.vjust=1,ticks= F))+##set colour
  theme(legend.title = element_text(colour = "black", size = rel(3.5)))+
  theme(strip.text = element_text(size = 45)) 

multiyear_Vmedium_plot <- add_coordinates(multiyear_Vmedium_plot, xmin = 70, xmax = 140, xstep = 15, ymin = 10, ymax = 60, ystep = 10, size = 1)

multiyear_Vmedium_plot_gt = ggplot_gtable(ggplot_build(multiyear_Vmedium_plot))
multiyear_Vmedium_plot_gt$heights = 1.3*multiyear_Vmedium_plot_gt$heights

CairoPNG(paste0("result/Fig/04.MISR_AOD_Volume/multiyear_volume_medium_plot.png"),width = 1800,height = 700, res=60)
print(grid.draw(multiyear_Vmedium_plot_gt))
dev.off()

##  plot multiyear volume ratio
colors <- RColorBrewer::brewer.pal(7,"RdYlGn")%>%rev()
breaks <- c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.001)

multiyear_V_ratio_result <- mutate(multiyear_V_ratio_result,
                                   discrete = cut(V_ratio, breaks = breaks, right = T, 
                                                         labels =  c("< 0.4","0.4 - 0.5","0.5 - 0.6","0.6 - 0.7","0.7 - 0.8", "0.8 - 0.9","> 0.9")))
summary(multiyear_V_ratio_result)
multiyear_V_ratio_plot <- grid_facet_plot(multiyear_V_ratio_result,  value = "discrete", facet = "year",facet_text = 18)+
  facet_wrap("year", nrow = 3, labeller = labeller(groupwrap = label_wrap_gen(10)))+
  geom_polygon(data=polygon,aes(x = long, y = lat, group = group),
               colour='black', fill='white',alpha=0.1)+
  scale_fill_manual(name = expression("Volume ratio"),
                    values =  colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2, 
                                         title.position="right",title.hjust=0.5,title.vjust=1,
                                         byrow = TRUE, ticks= F))+
  theme(legend.key.width=unit(1, "cm"))+
  theme(legend.key.height=unit(1, "cm"))+##legend height
  theme(legend.spacing.y = unit(.8, 'cm'))+
  theme(legend.title = element_text(colour = "black", size = rel(3.5)))+
  theme(strip.text = element_text(size = 35))

multiyear_V_ratio_plot <- add_coordinates(multiyear_V_ratio_plot, xmin = 70, xmax = 140, xstep = 15, ymin = 10, ymax = 60, ystep = 10, size = 1)

multiyear_V_ratio_plot_gt = ggplot_gtable(ggplot_build(multiyear_V_ratio_plot))
multiyear_V_ratio_plot_gt$heights = 1.3*multiyear_V_ratio_plot_gt$heights

CairoPNG(paste0("result/Fig/05.MISR_AOD_Volume_Ratio/multiyear_volume_ratio_plot.png"), width = 1800, height = 700, res = 60)
print(grid.draw(multiyear_V_ratio_plot_gt))
dev.off()

##-----06.Multiyear mean volume distribution plot---------
multiyear_mean_V_result <- fread(paste0(dir_predict,"Annual/MISR_GRID_particle_volume_multiyear_mean.csv"))

multiyear_mean_V_ratio_result <- multiyear_mean_V_result%>%
  group_by(x,y) %>% 
  summarise(V_ratio = V_mean[radius == 0.35]/V_mean[radius == 0.7])%>%
  filter(V_ratio != 0)
multiyear_mean_V_result
##--------06.1 Volume plot--------
colors <- RColorBrewer::brewer.pal(6,"YlGnBu")
multiyear_mean_V_plot  <- grid_facet_plot(data = multiyear_mean_V_result, value = "V_mean", facet = "radius", facet_text = 40)+
  facet_wrap("radius", nrow = 1, labeller = labeller(radius = c("0.35" = "Particle with radius less than 0.35 μm",
                                                                "0.7" = "Particle with radius less than 0.7 μm"),
                                                     groupwrap = label_wrap_gen(10)))+
  geom_polygon(data=polygon,aes(x = long, y = lat, group = group),
               colour='black', fill='white',alpha=0.1)+
  theme(legend.title = element_text(colour = "black", size = rel(2.5)))+
  scale_fill_gradientn(name=expression("Normalized Volume"),
                       colours = colors,
                       breaks=c(0,0.1,0.2,0.3,0.4,0.5),
                       labels=c(0,0.1,0.2,0.3,0.4,"> 0.5"),
                       trans = "sqrt",
                       guide = guide_colorbar(frame.colour = "black",frame.linewidth=2, 
                                              title.position="right",title.hjust=0.5, title.vjust=1,ticks= F)) ##set colour

multiyear_mean_V_plot <- add_coordinates(multiyear_mean_V_plot, xmin = 60, xmax = 150, xstep = 10, ymin = 0, ymax = 70, ystep = 10, size = 1.5)
multiyear_mean_V_plot
multiyear_mean_V_plot_gt = ggplot_gtable(ggplot_build(multiyear_mean_V_plot))
multiyear_mean_V_plot_gt$heights = 1.3*multiyear_mean_V_plot_gt$heights

CairoPNG(paste0("result/Fig/04.MISR_AOD_Volume/multiyear_mean_volume_plot.png"),width = 2800,height = 1500,res=110)
print(grid.draw(multiyear_mean_V_plot_gt))
dev.off()

##--------06.2 Ratio plot---------
colors <- RColorBrewer::brewer.pal(6,"RdYlGn")%>%rev()
breaks <- c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.001)
multiyear_mean_V_ratio_result <- mutate(multiyear_mean_V_ratio_result,
                                        discrete = cut(V_ratio, breaks = breaks, right = T,
                                                       labels =  c("< 0.4","0.4 - 0.5","0.5 - 0.6","0.6 - 0.7","0.7 - 0.8", "> 0.8")))

summary(multiyear_mean_V_ratio_result$V_ratio)
levels(multiyear_mean_V_ratio_result$discrete)
multiyear_mean_V_ratio_result
multiyear_mean_ratio_plot <- grid_plot(data = multiyear_mean_V_ratio_result, value = "discrete")+
  geom_polygon(data = polygon,aes(x = long, y = lat, group = group),
               colour='black', fill='white', alpha=0.1)+
  theme(legend.title = element_text(colour = "black", size = rel(2.5)))+
  scale_fill_manual(name = expression("Volume Ratio"),
                    values =  colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2, 
                                         title.position="right",title.hjust=0.5,title.vjust=1,
                                         byrow = TRUE, ticks= F))
multiyear_mean_ratio_plot <- add_coordinates(multiyear_mean_ratio_plot, xmin = 70, xmax = 140, xstep = 10, ymin = 10, ymax = 60, ystep = 10, size = 1.5)

multiyear_mean_ratio_plot

multiyear_mean_ratio_plot_gt = ggplot_gtable(ggplot_build(multiyear_mean_ratio_plot))
multiyear_mean_ratio_plot_gt$heights = 1.3*multiyear_mean_ratio_plot_gt$heights


CairoPNG(paste0("result/Fig/05.MISR_AOD_Volume_Ratio/multiyear_mean_volume_ratio_plot.png"),width = 2800,height = 1500,res=150)
print(grid.draw(multiyear_mean_ratio_plot_gt))
dev.off()

##-----07.Seasonal mean volume distribution plot---------
multiseason_mean_V_result <- readRDS(paste0(dir_predict, "Seasonal/MISR_GRID_particle_volume_multiseason_join.Rds"))
breaks <- c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.001)

multiseason_mean_V_ratio_result <- multiseason_mean_V_result%>%
  group_by(x, y, season, radius)%>%
  summarise(across("V_mean", mean))%>%
  summarise(V_ratio = V_mean[radius == 0.35]/V_mean[radius == 0.7])%>%
  mutate(discrete = cut(V_ratio, breaks = breaks, right = T, 
                        labels =  c("< 0.4","0.4 - 0.5","0.5 - 0.6","0.6 - 0.7","0.7 - 0.8","> 0.8")))
summary(multiseason_mean_V_ratio_result$V_ratio)
multiseason_mean_V_ratio_result
summary(multiseason_mean_V_ratio_result$V_ratio)
colors <- RColorBrewer::brewer.pal(6,"RdYlGn")%>%rev()
multiseason_ratio_plot  <- grid_facet_plot(data = multiseason_mean_V_ratio_result, value = "discrete", facet = "season", facet_text = 40)+
  facet_wrap("season", nrow = 2, labeller = labeller(groupwrap = label_wrap_gen(10)))+
  geom_polygon(data=polygon,aes(x = long, y = lat, group = group),
               colour='black', fill='white',alpha=0.1)+
  scale_fill_manual(name = expression("Volume Ratio"),
                    values =  colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2, 
                                         title.position="right",title.hjust=0.5,title.vjust=1,
                                         byrow = TRUE, ticks= F))+
  theme(legend.title = element_text(colour = "black", size = rel(2.5)))+
  theme(legend.key.width=unit(1, "cm"))+
  theme(legend.key.height=unit(1, "cm"))+##legend height
  theme(legend.spacing.y = unit(.8, 'cm'))+
  theme(strip.text = element_text(size = 35))

multiseason_ratio_plot <- add_coordinates(multiseason_ratio_plot, xmin = 60, xmax = 150, xstep = 10, ymin = 0, ymax = 70, ystep = 10, size = 1.5)
multiseason_ratio_plot_gt = ggplot_gtable(ggplot_build(multiseason_ratio_plot))
multiseason_ratio_plot_gt$heights = 1.3*multiseason_ratio_plot_gt$heights

CairoPNG(paste0("result/Fig/05.MISR_AOD_Volume_Ratio/multiseason_mean_volume_ratio_plot.png"),width = 2800,height = 2000,res=110)
print(grid.draw(multiseason_ratio_plot_gt))
dev.off()
