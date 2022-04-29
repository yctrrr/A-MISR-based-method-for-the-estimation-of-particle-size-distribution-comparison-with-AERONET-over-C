library(sp)
library(dplyr)
library(rgdal)
library(data.table)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(Cairo)
library(stringr)

source("Point_visulaization.R")
##------01.Plot Aeronet AOD station----------
##  China polygon
polygon <- readOGR(dsn="data/other_data/prov_China/bou2_4p_region.shp",
                   layer="bou2_4p_region")
polygon@proj4string <-  CRS("+proj=longlat +datum=WGS84")
polygon_sf <- st_as_sf(polygon)

##  read in AERONET data during 2004-2016
aeronet <- fread(paste0("data/AOD/Aeronet AOD/aeronet_aod_2004-2016.out"))%>%
  rename(Lon="Longitude(Degrees)",Lat="Latitude(Degrees)", AERONET_Site="AERONET_AERONET_Site")%>%
  mutate(AERONET_Site = str_replace_all(AERONET_Site, "_|-", " "))

aeronet_sp <- SpatialPointsDataFrame(coords = aeronet[,c("Lon","Lat")],
                                     data = aeronet,
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))

aeronet_china <- aeronet[!is.na(sp::over(aeronet_sp,polygon)[[1]]),]

## AERONET sites names
aeronet_china$AERONET_Site%>%unique()
aeronet_china_match <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_annual_distribution.Rds"))%>%
  count(AERONET_Site,Lon,Lat)
sites_all <- unique(aeronet_china_match$AERONET_Site)

##  label main sites
monthly_count <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_monthly_distribution.Rds"))%>%
  count(AERONET_Site,Lon,Lat)
sites_sel <- monthly_count[monthly_count$n >= 792,"AERONET_Site"]

aeronet_china_count <- aeronet_china%>%
  count(AERONET_Site, Lon, Lat)%>%
  filter(AERONET_Site%in%sites)

aeronet_china_label <- filter(aeronet_china_count, AERONET_Site%in%sites_sel)

# distinct(aeronet_china,AERONET_AERONET_Site, Lon, Lat, .keep_all = TRUE)%>%nrow()
aeronet_china_sp <- SpatialPointsDataFrame(coords = aeronet_china_count[,c("Lon","Lat")],
                                           data = aeronet_china_count,
                                           proj4string = CRS("+proj=longlat +datum=WGS84"))

aeronet_china_sf <- st_as_sf(aeronet_china_sp)%>%
  dplyr::rename(Daily_observation_count = n)


##  write in AERONET sites information
fwrite(aeronet_china_sf, file = "result/Table/Summary/Aeronet_main_sites.csv")
aeronet_china_sf
##  visualize AERONET sites
colors <- brewer.pal(6, "GnBu")
aeronet_china_sf$Daily_observation_count%>%summary()
aeronet_site_plot <- sf_point_plot(aeronet_china_sf, "Daily_observation_count",polygon_sf)+
  scale_fill_gradientn(name = expression("Daily observation count"),
                       colours = colors, trans = "log10",
                       limits = c(1,2000),
                       guide = guide_colorbar(frame.colour = "black", frame.linewidth = 1,
                       title.position = "right", title.hjust = 0.5, title.vjust = 1, ticks = F))+
  geom_label_repel(data = aeronet_china_label ,
                  aes(x = Lon, y = Lat, label = AERONET_Site),
                  label.size = NA, alpha = 1, label.padding=.1,
                  point.padding = unit(0.5, "lines"),
                  na.rm=TRUE, seed = 1234, 
                  size = 5, box.padding = unit(0.6, "lines") # size and size distance
                  ) ## add aeronet sites

aeronet_site_plot
ggsave("result\\Fig\\01.Data_analysis\\Aeronet_site_distribution.png",aeronet_site_plot,width = 12,height = 7)

#-------02.Radius and distribution-----------
## Radius range of 8 MISR AOD components
range1 <- c(0.001,0.4)
range2 <- c(0.001,0.75)
range3 <- c(0.01,1.5)
range6 <- c(0.1,50)
range8 <- c(0.001,0.75)
range14 <- c(0.001,0.75)
range19 <- c(0.1,1)
range21 <- c(0.1,6)
##  Distribution width
Dw <- c(1.65,1.7,1.75,1.9,1.7,1.7,1.5,2)
##  Volume-weighted mode radius of 8 components
# radius_list <- c(0.05,0.105,0.225,2.275,0.105,0.105,0.7,2.61)

##  Number-weighted mode radius of 8 components
radius_list <- c(0.03,0.06,0.12,1,0.06,0.06,0.5,1)

##  sequence number of MISR AOD components
number <- c(1,2,3,6,8,14,19,21)

##  function to reconstruct data frame
framer <-function(data){return  (data.frame(y=data$y, x=data$x)) }

##  Main input: mean radius and distribution width
for (i in 1:length(number)) {
  ##  number-weighted distribution
  ND <- rlnorm(50000, meanlog = log(radius_list[i]), sdlog = log(Dw[i]))

  ##  calculate density of lognormal distribution
  range <- get(paste0("range", number[i]))
  density <- density(ND, from=0, to=50, n=5000)#n=5000

  ##  convert to datafrmae
  dens_df <- framer(density)

  ##  density outside radius range is 0
  sel <- dens_df$x<range[1]|dens_df$x>range[2]
  dens_df$y[sel] <- 0

  ##  rename MISR component composition
  if(i==2|i==5|i==6){
    dens_df$component <- paste0("Component 2/8/14")
  }else{
    dens_df$component <- paste0("Component ",number[i])
  }
  ##  asign value name
  assign(paste0("d",number[i]),dens_df)
}
save(list = paste0("d",number),file = "result/Rdata/MISR_AOD_Parameter/MISR_AOD_SD.Rdata")

dens_df_full <- do.call(rbind,mget(c("d1","d2","d3","d6","d19","d21")))
dens_df_full$component<- factor(dens_df_full$component,levels=paste0("Component ",c("1","2/8/14","3","6","19","21")))
dens_df_full
summary(dens_df$x)
##  GGplot
##  mapping colour
color <-  RColorBrewer::brewer.pal(6, 'Set1')
color[6] <- "#d8b365"
component_SD <- ggplot()+
  geom_line(data=dens_df_full,aes(x = x,y = y,color = component),size=0.8)+
  scale_x_continuous(limits = c(0, 1.5),breaks = seq(0, 1.5,0.2))+
  labs(x= expression("Number Weighted Particle Radius(µm)"),
       y= "Mass Fraction(dN/dr)",color=NULL,vjust=1)+
  scale_color_manual(values =  color)+
  theme_bw()+
  theme(legend.position=c(0.8,0.6))+
  theme(legend.key.width=unit(2.5, "cm"))+
  theme(legend.key.height=unit(2, "cm"))+##legend height
  theme(legend.text=element_text(colour = "black", size = rel(3)))+ ##legend text size
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=40,vjust=1.3))+
  theme(panel.border = element_rect(colour = "black",size=1))##panel line width

CairoPNG("result/Fig/01.Data_analysis/MISR_radius_plot.png",width = 1500,height = 1000,res = 70)
component_SD
dev.off()

