library(sp)
library(raster)
library(dplyr)
library(rgdal)
library(data.table)
library(sf)
library(stringr)
library(geosphere)

source("Grid_visualization.R")
source("Spatial_analysis.R")

##  China polygon
polygon <- readOGR(dsn="data/other_data/prov_China/bou2_4p_region.shp",
                   layer="bou2_4p_region")
polygon@proj4string <-  CRS("+proj=longlat +datum=WGS84")

year <- 2016
##------01.Daily matching------
for (year in 2004:2016) {
  ##  read in observed MISR AOD data
  aod <- fread(paste0("data/AOD/MISR AOD/China_MISR_AOT_component_AGP_",year,".csv"))
  ##  Calculate total AOD
  aod <- mutate(aod, Component_all = Component_01 + Component_02 + Component_03 + Component_06 +
                  Component_08 + Component_14 + Component_19 + Component_21)

  ##  read in aeronet data
  aeronet <- fread(paste0("data/AOD/Aeronet AOD/aeronet_aod_",year,".out"))%>%
    rename(Lon="Longitude(Degrees)",Lat="Latitude(Degrees)")
  
  aeronet_sp <- SpatialPointsDataFrame(coords = aeronet[,c("Lon","Lat")],
                                       data = aeronet,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

  T.lst<- c(seq(as.Date(paste(year,"/1/1",sep='')),
                as.Date(paste(year,"/12/31",sep='')), by = "day"))
  rm(merge_data)
  
  flag = 0
  for (i in 1:length(T.lst)) {
    print(i)
    ##  read in daily MISR AOD data
    aod_sel <- filter(aod,DOY==i)%>%
      filter(!(Component_01==0&Component_02==0&Component_03==0&Component_06==0&
                 Component_08==0&Component_14==0&Component_19==0&Component_21==0))
    

    if(nrow(aod_sel)>0){
      aod_sp <- SpatialPointsDataFrame(coords = aod_sel[,c("MISRLon","MISRLat")],
                                       data = aod_sel,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

      ##  extract aeronet AOD in china
      aeronet_china <- aeronet[!is.na(over(aeronet_sp,polygon)[[1]]),]
      aeronet_china_sel <- aeronet_china[grepl(year, aeronet_china$`Date(dd:mm:yyyy)`, fixed=TRUE),]%>%
        filter(Day_of_Year==i)

      # unique(aeronet_china_sel$AERONET_AERONET_Site)
      if(nrow(aeronet_china_sel) > 0){
        # plot(aeronet_china_sp)
        aeronet_china_site <- aeronet_china_sel%>%
          dplyr::select(c(1,6:27,"Lon","Lat"))%>%
          group_by(AERONET_AERONET_Site, Lon, Lat)%>%
          summarise_all(mean)%>%
          as.data.frame()
    
        aeronet_china_sp <- SpatialPointsDataFrame(coords = aeronet_china_site[,c("Lon","Lat")],
                                                   data = aeronet_china_site,
                                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
        ##  buffer within 10km
        merge_data <- circle_buffer(data = aeronet_china_sp, match_data = aod_sp[,c(15:16,18,22:29)])

        if(nrow(merge_data)>0) {
          if(flag == 0){
            valid_temp <- merge_data
            flag = 1
            print("success")
          }else if(flag == 1){
            valid_temp <- full_join(valid_temp, merge_data)
          }
        }
      }
    }
  }## end of day
  valid_df_daily <- valid_temp
  saveRDS(valid_df_daily,file = paste0("result/Rdata/Aerronet_AOD_valid/Daily/aerronet_vs_misr_daily_",year,".Rds"))
}
year = 2016
readRDS(paste0("result/Rdata/Aerronet_AOD_valid/Daily/aerronet_vs_misr_daily_",year,".Rds"))

##------02.Monthly matching------
for (year in 2004:2016) {
  ##  read in observed MISR AOD data
  aod <- fread(paste0("data/AOD/MISR AOD/China_MISR_AOT_component_AGP_",year,".csv"))
  ##  Calculate total AOD
  aod <- mutate(aod, Component_all = Component_01 + Component_02 + Component_03 + Component_06 +
                  Component_08 + Component_14 + Component_19 + Component_21)
  ##  read in aeronet data
  aeronet <- fread(paste0("AOD/Aeronet AOD/aeronet_aod_",year,".out"))%>%
    rename(Lon="Longitude(Degrees)",Lat="Latitude(Degrees)")
  aeronet_sp <- SpatialPointsDataFrame(coords = aeronet[,c("Lon","Lat")],
                                       data = aeronet,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))
  T.lst<- c(seq(as.Date(paste(year,"/1/1",sep='')),
                as.Date(paste(year,"/12/31",sep='')), by = "day"))
  
  flag = 0
  rm(valid_df_monthly)
  for (i in 1:12) {
    print(i)
    ##  read in daily MISR AOD data
    aod_sel <- filter(aod,Month==i)%>%
      filter(!(Component_01==0&Component_02==0&Component_03==0&Component_06==0&
                 Component_08==0&Component_14==0&Component_19==0&Component_21==0))%>%
      filter(AOT_bestestimate!= -9999)
  
    if(nrow(aod_sel)>0){
      aod_sp <- SpatialPointsDataFrame(coords = aod_sel[,c("MISRLon","MISRLat")],
                                       data = aod_sel,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))
  
      ##  extract aeronet AOD in china
      aeronet_china <- aeronet[!is.na(over(aeronet_sp,polygon)[[1]]),]
      aeronet_china$Year <- substr(aeronet_china$`Date(dd:mm:yyyy)`,7,10)
      aeronet_china$date <- as.Date(aeronet_china$Day_of_Year-1, origin = paste0(aeronet_china$Year,"-01-01"))
      aeronet_china$month <- month(aeronet_china$date)
      aeronet_china_sel <- filter(aeronet_china,Year==year&month==i)
  
      if(nrow(aeronet_china_sel)>0){
        aeronet_china_site <- aeronet_china_sel%>%
          dplyr::select(c(1,6:27,"Lon","Lat"))%>%
          group_by(AERONET_AERONET_Site, Lon, Lat)%>%
          summarise_all(mean)%>%
          as.data.frame()
  
        aeronet_china_sp <- SpatialPointsDataFrame(coords = aeronet_china_site[,c("Lon","Lat")],
                                                   data = aeronet_china_site,
                                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
        
        merge_data <- circle_buffer(data = aeronet_china_sp, match_data = aod_sp[,c(15:16,22:29)])
        
        if(nrow(merge_data)>0) {
          if(flag == 0){
            valid_temp <- merge_data
            flag = 1
            print("success")
          }else if(flag == 1){
            valid_temp <- full_join(valid_temp, merge_data)
          }
        }
      }
    }
  }
  valid_df_monthly <- valid_temp
  saveRDS(valid_df_monthly,file = paste0("result/Rdata/Aerronet_AOD_valid/Monthly/aerronet_vs_misr_monthly_",year,".Rds"))
}

##------03.Seasonal matching------
for (year in 2004:2016) {
  season_name <- c("Spring","Summer","Autumn","Winter")
  ##  read in observed MISR AOD data
  aod <- fread(paste0("data/AOD/MISR AOD/China_MISR_AOT_component_AGP_",year,".csv"))
  
  ##  Calculate total AOD and identify season ids
  aod <- mutate(aod, Component_all = Component_01 + Component_02 + Component_03 + Component_06 +
                  Component_08 + Component_14 + Component_19 + Component_21)%>%
    mutate(season = ifelse(Month%in%c(3,4,5),"Spring", 
                  ifelse(Month%in%c(6,7,8),"Summer",
                         ifelse(Month%in%c(9,10,11),"Autumn","Winter"))))%>%
    mutate(across(season,factor,levels = season_name))

  aod$season_id <- match(aod$season,season_name)
  
  ##  read in aeronet data
  aeronet <- fread(paste0("data/AOD/Aeronet AOD/aeronet_aod_",year,".out"))%>%
    rename(Lon="Longitude(Degrees)",Lat="Latitude(Degrees)")

  aeronet_sp <- SpatialPointsDataFrame(coords = aeronet[,c("Lon","Lat")],
                                       data = aeronet,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))
  T.lst<- c(seq(as.Date(paste(year,"/1/1",sep='')),
                as.Date(paste(year,"/12/31",sep='')), by = "day"))
  
  flag = 0
  rm(valid_df_seasonal)
  for (i in 1:4) {
    print(i)
    ##  read in daily MISR AOD data
    aod_sel <- filter(aod,season_id==i)%>%
      filter(!(Component_01==0&Component_02==0&Component_03==0&Component_06==0&
                 Component_08==0&Component_14==0&Component_19==0&Component_21==0))%>%
      filter(AOT_bestestimate!= -9999)

    if(nrow(aod_sel)>0){
      aod_sp <- SpatialPointsDataFrame(coords = aod_sel[,c("MISRLon","MISRLat")],
                                       data = aod_sel,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

      ##  extract aeronet AOD in china
      aeronet_china <- aeronet[!is.na(over(aeronet_sp,polygon)[[1]]),]
      aeronet_china$Year <- substr(aeronet_china$`Date(dd:mm:yyyy)`,7,10)
      aeronet_china$date <- as.Date(aeronet_china$Day_of_Year-1, origin = paste0(aeronet_china$Year,"-01-01"))
      aeronet_china$month <- month(aeronet_china$date)

      ##  seasonal aeronet aod data
      aeronet_china <- aeronet_china%>%
        mutate(season = as.factor(ifelse(month%in%c(3,4,5),"Spring", 
                               ifelse(month%in%c(6,7,8),"Summer",
                                      ifelse(month%in%c(9,10,11),"Autumn","Winter")
                                      ))))%>%
        mutate(across(season,factor,levels = season_name))
      
      aeronet_china$season_id <- match(aeronet_china$season,season_name)
      aeronet_china_sel <- filter(aeronet_china,Year==year&season_id==i)
      
      if(nrow(aeronet_china_sel)>0){
        aeronet_china_site <- aeronet_china_sel%>%
          dplyr::select(c(1,6:27,"Lon","Lat"))%>%
          group_by(AERONET_AERONET_Site, Lon, Lat)%>%
          summarise_all(mean)%>%
          as.data.frame()
      
        aeronet_china_sp <- SpatialPointsDataFrame(coords = aeronet_china_site[,c("Lon","Lat")],
                                                   data = aeronet_china_site,
                                                   proj4string = CRS("+proj=longlat +datum=WGS84"))

        ##  buffer within 10km
        merge_data <- circle_buffer(data = aeronet_china_sp, match_data = aod_sp[,c(15:16,18,22:29)])

        if(nrow(merge_data)>0) {
          if(flag == 0){
            valid_temp <- merge_data
            flag = 1
            print("success")
          }else if(flag == 1){
            valid_temp <- full_join(valid_temp, merge_data)
          }
        }
      }
    }
  }
  valid_df_seasonal <- valid_temp
  saveRDS(valid_df_seasonal,file = paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_seasonal_",year,".Rds"))
}

valid_df_seasonal <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_seasonal_",year,".Rds"))
valid_df_seasonal

##------04.Annual matching------
for (year in 2004:2016) {
  ##  read in observed MISR AOD data
  aod <- fread(paste0("data/AOD/MISR AOD/China_MISR_AOT_component_AGP_",year,".csv"))%>%
    mutate(Component_all = Component_01 + Component_02 + Component_03 + Component_06 +
                  Component_08 + Component_14 + Component_19 + Component_21)
  ##  read in aeronet data
  aeronet <- fread(paste0("D:/shaoyanchuan/project/pm_distribution/data/AOD/Aeronet AOD/aeronet_aod_",year,".out"))%>%
    rename(Lon="Longitude(Degrees)",Lat="Latitude(Degrees)")
  
  aeronet_sp <- SpatialPointsDataFrame(coords = aeronet[,c("Lon","Lat")],
                                       data = aeronet,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))
  flag = 0
  rm(valid_df_annual)
  ##  read in daily MISR AOD data
  aod_sel <- filter(aod)%>%
    filter(!(Component_01==0&Component_02==0&Component_03==0&Component_06==0&
               Component_08==0&Component_14==0&Component_19==0&Component_21==0))%>%
    filter(AOT_bestestimate!= -9999)

  aod_sp <- SpatialPointsDataFrame(coords = aod_sel[,c("MISRLon","MISRLat")],
                                   data = aod_sel,
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))

  ##  extract aeronet AOD in china
  aeronet_china <- aeronet[!is.na(over(aeronet_sp,polygon)[[1]]),]
  aeronet_china$Year <- substr(aeronet_china$`Date(dd:mm:yyyy)`,7,10)
  aeronet_china_sel <- filter(aeronet_china,Year==year)
  
  aeronet_china_site <- aeronet_china_sel%>%
    dplyr::select(c(1,6:27,"Lon","Lat"))%>%
    group_by(AERONET_AERONET_Site, Lon, Lat)%>%
    summarise_all(mean)%>%
    as.data.frame()
      
  aeronet_china_sp <- SpatialPointsDataFrame(coords = aeronet_china_site[,c("Lon","Lat")],
                                             data = aeronet_china_site,
                                             proj4string = CRS("+proj=longlat +datum=WGS84"))
  ##  buffer within 10km
  merge_data <- circle_buffer(data = aeronet_china_sp, match_data = aod_sp[,c(15,22:29)])
   
  if(nrow(merge_data)>0) {
     if(flag == 0){
       valid_temp <- merge_data
       flag = 1
       print("success")
     }else if(flag == 1){
      valid_temp <- full_join(valid_temp, merge_data)
   }
  }
  saveRDS(valid_df_annual, file = paste0("result/Rdata/Aerronet_AOD_valid/Annual/aerronet_vs_misr_annual_",year,".Rds"))
}
