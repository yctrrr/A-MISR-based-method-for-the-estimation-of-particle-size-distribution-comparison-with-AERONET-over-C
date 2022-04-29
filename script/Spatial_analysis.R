library(sf)
library(dplyr)
library(sp)
library(data.table)
library(raster)
##  normalize function
standard <- function(x){(x-min(x))/(max(x)-min(x))}


##-----------Convert spatial points to grids: mean value---------------
points_to_grids <- function(data, x = "x", y = "y", crs = "+proj=longlat +datum=WGS84", polygon, raster){
  start_time <- Sys.time()
  # data_sp <- SpatialPointsDataFrame(select(data, x, y), data, proj4string = CRS(crs))
  # data_poly <- data_sp[!is.na(over(data_sp,polygon)[[1]]),]
  # data_sf_poly <- st_as_sf(data_poly, coords = c(x, y), crs = crs)
  # data_raster_poly <- rasterize(data_sf_poly, raster, fun = mean)

  data_sf <- st_as_sf(data, coords = c(x, y), crs = crs)
  data_raster <- rasterize(data_sf, raster , fun = mean)
  data_raster_poly <- mask(data_raster, polygon)
  # endtime <- Sys.time()
  # print(endtime - start_time)
  data_griddf_poly <- data_raster_poly%>%
    rasterToPoints()%>%
    as.data.frame()
  return(data_griddf_poly)
}


##--------mask points outside polygons----------
points_mask <- function(data, x = "x", y = "y", crs = "+proj=longlat +datum=WGS84", polygon){
  data_sp <- SpatialPointsDataFrame(dplyr::select(data, x, y), data, proj4string = CRS(crs))
  data_sp_poly <- data_sp[!is.na(over(data_sp,polygon)[[1]]),]
  return(data_sp_poly@data)
}

##-----------Circle buffer average value calculation---------------
##  data, match_data: spatial points data
##  buffer: searching radius, default 10km
## ...: add match variable names
circle_buffer <- function(data, match_data, buffer = 10000){
  for (p in 1:nrow(data)) {
    data_sel <- data[p,]
    match_data_all <- match_data[distHaversine(data_sel, match_data) < buffer,]@data
    if(nrow(match_data_all) > 0){
      match_data_mean <- colMeans(match_data_all, na.rm = T)%>%t()%>%as.data.frame()
    }else{
      match_data_mean <- match_data_all[NA,]
    }
    merge_data_temp <- cbind(data_sel@data, match_data_mean)
    if(p == 1){
      merge_data <- merge_data_temp
    }else{
      merge_data <- full_join(merge_data,merge_data_temp)
    }
  }
  return(merge_data)
}
