library(sp)
library(raster)
library(dplyr)
library(rgdal)
library(data.table)
library(sf)
library(stringr)
library(geosphere)
library(ggh4x)
library(ggplot2)

source("Grid_visualization.R")
source("Spatial_analysis.R")
source("Scatter_plot.R")
source("AOD_to_VSD.R")

##*  01.Sample-based VSD calculation--------
match_level <- c("daily","monthly","seasonal","annual")
##  normalized function
standard <- function(x){(x-min(x))/(max(x)-min(x))}
level = 1
##  calculate MISR VSD and match them with AERONET
for (level in 1:length(match_level)) {
  flag = 0
  for (year in 2004:2016) {
    valid_df <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/",str_to_title(match_level[level]),"/aerronet_vs_misr_",match_level[level],"_",year,".Rds"))
    bin_sel <- c(4:25)

    rm_valid <- colnames(valid_df)[bin_sel]%>%
      gsub("X","",.)%>%
      as.numeric()

    valid_df_obs <- reshape(valid_df, varying = bin_sel, v.names = "VD_ln_obs", 
            timevar = "radius", times = rm_valid, direction = "long")%>%
      mutate(VD_obs = VD_ln_obs/radius)
    rownames(valid_df_obs) <- NULL
    
    valid_df_calcu <- misr_aod_to_distri(valid_df, rm_valid = rm_valid)%>%
      filter(radius%in%rm_valid)
    
    valid_df_reshape <- left_join(valid_df_calcu, valid_df_obs)%>%na.omit()

    cor(valid_df_reshape$VD, valid_df_reshape$VD_obs)
    valid_ID <- unique(valid_df_reshape$AERONET_AERONET_Site)

    ##  combine
    if(flag == 0){
      distri_valid_allyear <- valid_df_reshape
      flag = 1
    }else{
      distri_valid_allyear <- full_join(distri_valid_allyear,valid_df_reshape)
    }
  }
  distri_valid_allyear <- mutate(distri_valid_allyear, VD = standard(distri_valid_allyear$VD), VD_obs = standard(distri_valid_allyear$VD_obs))%>%
    dplyr::rename(AERONET_Site = AERONET_AERONET_Site)%>%
    mutate(AERONET_Site = str_replace_all(AERONET_Site, "_|-", " "))
  saveRDS(distri_valid_allyear, paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_",match_level[level],"_distribution.Rds"))
}


##---------*02. Scatter plot of dV/dr for different scales----------
library(Cairo)
match_level <- c("daily","monthly","seasonal","annual")
count <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_monthly_distribution.Rds"))%>%
  count(AERONET_Site,Lon,Lat)
count
sites_sel <- count[count$n >= 792,"AERONET_Site"]
sites_sel
level = 2
for (level in 1:length(match_level)) {
  distri_valid_allyear <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_",match_level[level],"_distribution.Rds"))%>%
    mutate(Site_adjust = ifelse(AERONET_Site%in%sites_sel,AERONET_Site,"Other sites"),
           # discrete = cut(radius, breaks = c(0, 0.1, 0.2, 0.35, 0.7, 3, +Inf), labels = c("< 0.1","0.1 - 0.2","0.2 - 0.35", "0.35 - 0.7", "0.7 - 3", "> 0.3")))
           discrete = cut(radius, breaks = c(0, 0.1, 0.3, 0.8, 2.5, 8, +Inf), labels = c("< 0.1","0.1 - 0.3","0.3 - 0.8", "0.8 - 2.5", "2.5 - 8", "> 8")))
           
  distri_valid_allyear$AERONET_Site%>%unique()
  distri_valid_radius <- distri_valid_allyear%>%
    full_join(mutate(filter(distri_valid_allyear, radius <= 0.7), discrete = "< 0.7"))%>%
    full_join(mutate(filter(distri_valid_allyear, radius <= 0.35), discrete = "< 0.35"))%>%
    mutate(discrete = factor(discrete, levels = c("< 0.1","0.1 - 0.3","0.3 - 0.8", "0.8 - 2.5", "2.5 - 8", "> 8","< 0.35","< 0.7")))

  distri_valid_radius$discrete
  site_name <- unique(distri_valid_allyear$Site_adjust)
  site_order <- sort(site_name[!site_name%in%"Other sites"])%>%
    c("Other sites")
  
  allyear_df <- data.frame(y = distri_valid_allyear$VD_obs, ypred = distri_valid_allyear$VD, Year = distri_valid_allyear$Year,
                           Site = distri_valid_allyear$Site_adjust, radius_label = distri_valid_allyear$discrete)%>%
    mutate(across(Site, factor, levels=site_order))

  allyear_df_radius <- data.frame(y = distri_valid_radius$VD_obs, ypred = distri_valid_radius$VD, Year = distri_valid_radius$Year,
                           Site = distri_valid_radius$Site_adjust, radius_label = distri_valid_radius$discrete)%>%
    mutate(across(Site, factor, levels=site_order))
  
  ##---------**Scatter plot of dV/dr for overall----------
  xlim = ylim = c(0, 1)
  allyear_validation <- cv.plot(allyear_df, xlim =xlim,ylim=ylim,
          cex = 10, alpha = 0.3, binwidth=c(1.2, 1.2),ylimrange = 0.075)+
    stat_binhex(bins = 150)+
    scale_fill_gradientn(name = expression("Count"),
                         trans = "log10",
                         colours = c("#114bf9","#0ffff4","#60f818",
                                     "#f2f607","#f2d01e","#f20505"),
                         guide = guide_colorbar(frame.colour = "black",frame.linewidth=2,
                                                title.position="top",title.hjust=0.5,title.vjust = 1,ticks= F))+
    labs(x = expression(atop("AERONET derived dV/dr")),
         y = expression("MISR derived dV/dr"))+
    theme(legend.position="right")+
    theme(legend.key.width=unit(.6, "cm"))+
    theme(legend.key.height=unit(2.8, "cm"))+##legend height
    theme(legend.text=element_text(colour = "black", size = rel(2)))+ ##legend text size
    theme(legend.title = element_text(colour = "black", size = rel(2.5)))+
    theme(legend.margin=margin(r= 10))+
    theme(plot.subtitle = element_text(colour = "black", size = rel(2.8),hjust = 0.5))+
    theme(axis.text=element_text(size=24), axis.title=element_text(size=25))
  
  CairoPNG(paste0("result/Fig/02.AOD_scatter_validation/allyear_",match_level[level],"_validation.png"),width = 3200,height = 2800,res=300)
  print(allyear_validation)
  dev.off()

  ##---------**Scatter plot of dV/dr for each site----------
  ##  multiple plots for each site
  # xlim = ylim = NULL
  xlim = ylim = c(0, 1)
  allyear_validation_facet <- cv.plot.facet.blue(allyear_df, nrow = 2, xlim = xlim,ylim = ylim,
                                                 cex = 6, alpha = 0.3, coord_ylim=1,
                                                 ylimrange = 0.0015,facet = "Site",label = "", vjust = 1.1)+
    labs(x = expression(atop("AERONET derived dV/dr")),
         y = expression("MISR derived dV/dr"))+
    theme(legend.position = "right")+
    theme(legend.key.width = unit(.6, "cm"))+
    theme(legend.key.height = unit(2.8, "cm"))+##legend height
    theme(legend.text = element_text(colour = "black", size = rel(3)))+ ##legend text size
    theme(legend.title = element_text(colour = "black", size = rel(3.5)))+
    theme(legend.margin=margin(r= 10))+
    theme(plot.subtitle = element_text(colour = "black", size = rel(5),hjust = 0.5))+
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 30))+
    theme(panel.spacing.x = unit(2, "lines"))
  
  CairoPNG(paste0("result/Fig/02.AOD_scatter_validation/allyear_",match_level[level],"_validation_each_site.png"),width = 4200,height = 2500,res=280)
  print(allyear_validation_facet)
  dev.off()

  ##---------**Scatter plot of dV/dr for each radius range----------
  xlim = ylim = c(0, 1)
  radius_validation <- cv.plot.facet.blue(allyear_df_radius, nrow = 2,xlim = xlim,ylim = ylim, 
                                          cex = 7, alpha = 0.3, binwidth=c(1.2, 1.2),ylimrange = 0.075,
                                          facet = "radius_label",label = "Radius ",vjust = 1.1)+
    stat_binhex(bins = 150)+
    scale_fill_gradientn(name = expression("Count"),
                         trans = "log10",
                         colours = c("#114bf9","#0ffff4","#60f818",
                                     "#f2f607","#f2d01e","#f20505"),
                         guide = guide_colorbar(frame.colour = "black",frame.linewidth=2,
                                                title.position="top",title.hjust=0.5,title.vjust = 1,ticks= F))+
    labs(x = expression(atop("AERONET derived dV/dr")),
         y = expression("MISR derived dV/dr"))+
    theme(legend.position="right")+
    theme(legend.key.width=unit(.6, "cm"))+
    theme(legend.key.height=unit(2.8, "cm"))+##legend height
    theme(legend.text=element_text(colour = "black", size = rel(2)))+ ##legend text size
    theme(legend.title = element_text(colour = "black", size = rel(2.5)))+
    theme(legend.margin=margin(r= 10))+
    theme(plot.subtitle = element_text(colour = "black", size = rel(2.5),hjust = 0.5))+
    theme(axis.text=element_text(size=20), axis.title=element_text(size=30)) +
    theme(panel.spacing.x = unit(2, "lines"))

  CairoPNG(paste0("result/Fig/02.AOD_scatter_validation/radius_range_",match_level[level],"_validation.png"),width = 5000,height = 2500,res=240)
  print(radius_validation)
  dev.off()
  
  ##---------**Scatter plot of dV/dr for each year----------
  eachyear_validation <- cv.plot.facet.blue(allyear_df,nrow = 3,xlim = xlim,ylim = ylim, 
                                            cex = 7, alpha = 0.3, binwidth=c(1.2, 1.2),ylimrange = 0.075,vjust = 1.1)+
    labs(x = expression(atop("AERONET derived dV/dr")),
         y = expression("MISR derived dV/dr"))+
    theme(legend.position = "right")+
    theme(legend.key.width = unit(.6, "cm"))+
    theme(legend.key.height = unit(2.8, "cm"))+##legend height
    theme(legend.text = element_text(colour = "black", size = rel(3)))+ ##legend text size
    theme(legend.title = element_text(colour = "black", size = rel(3.5)))+
    theme(legend.margin=margin(r= 10))+
    theme(plot.subtitle = element_text(colour = "black", size = rel(5),hjust = 0.5))+
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 35))+
    theme(panel.spacing.x = unit(2, "lines"))

  CairoPNG(paste0("result/Fig/02.AOD_scatter_validation/eachyear_",match_level[level],"_validation.png"),width = 5800,height = 3800,res=280)
  print(eachyear_validation)
  dev.off()
}

##---------* 03.Distribution validation----------
valid_df <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_annual_distribution.Rds"))
bin_sel <- c(4:25)
rm_valid <- colnames(valid_df)[bin_sel]%>%
  gsub("X","",.)%>%
  as.numeric()
count <- count(valid_df, AERONET_Site,Lon,Lat)
site_name <- count[count$n >= 660,"AERONET_Site"]%>%
  c("Other sites")
site_order <- sort(site_name[!site_name%in%"Other sites"])%>%
  c("Other sites")

label <- seq(0, 15, by=15/21)
plot_label <- formatC(rm_valid[seq(1,22,3)], width = 2, digits = 2, format = "f")
breaks <- label[seq(1,22,3)]

distri_valid_allyear <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_daily_distribution.Rds"))%>%
  mutate(Site_adjust = ifelse(AERONET_Site%in%sites_sel,AERONET_Site,"Other sites"))%>%
  mutate(across(Site_adjust, factor, levels = site_order))%>%
  rename(MISR = VD)%>%
  rename(AERONET = VD_obs)

distri_valid_mean <- distri_valid_allyear%>%
  group_by(radius)%>%
  summarise(across(c("MISR","AERONET"), mean))%>%
  as.data.frame()%>%
  melt(measure.vars=c("MISR","AERONET"))
label <- seq(0,15, by=15/21)
distri_valid_mean$label <- label

distri_valid_each_mean <- distri_valid_allyear%>%
  group_by(radius, Year)%>%
  summarise(across(c("MISR","AERONET"), mean))%>%
  as.data.frame()%>%
  melt(measure.vars=c("MISR","AERONET"))
label <- rep(seq(0, 15, by=15/21), each = 13)
distri_valid_each_mean$label <- label

distri_compare <- curve.plot.facet(distri_valid_mean)+
  labs(x = "Radius(μm)",y = "Derived dV/dr")+
  scale_x_continuous(breaks = breaks, labels = plot_label)+
  theme(legend.key.width=unit(1, "cm"))+
  theme(legend.key.height=unit(1, "cm"))+##legend height
  theme(legend.spacing.y = unit(.8, 'cm'))+
  theme(legend.text=element_text(colour = "black", size = rel(1)))+ ##legend text size
  theme(legend.title = element_blank())+
  theme(strip.text = element_text(size = rel(1.2)))+
  theme(plot.subtitle = element_text(colour = "black", size = rel(1.4),hjust = 0.5))+
  theme(axis.text=element_text(size=rel(1.2)),
        axis.title=element_text(size=rel(1.2)))


CairoPNG(paste0("result/Fig/03.AOD_distribution_validation/aeronet_vs_misr_multiyear_mean.png"),
         width = 3000,height = 1400,res=300)
print(distri_compare)
dev.off()

distri_compare_each <- curve.plot.facet(distri_valid_each_mean)+
  facet_nested_wrap( ~  Year + variable, ncol = 10, scales = "free")+
  labs(x = "Radius",y = "Derived dV/dr")+
  theme(legend.key.width=unit(2, "cm"))+
  theme(legend.key.height=unit(2, "cm"))+##legend height
  theme(legend.spacing.y = unit(2, 'cm'))+
  theme(legend.text=element_text(colour = "black", size = rel(1.8)))+ ##legend text size
  theme(legend.title = element_blank())+
  theme(strip.text = element_text(size = rel(1.8)))+
  theme(plot.subtitle = element_text(colour = "black", size = rel(2),hjust = 0.5))+
  theme(axis.text=element_text(size=rel(2)),
        axis.title=element_text(size=rel(2)),axis.title.x = element_text(vjust = -0.1))

CairoPNG(paste0("result/Fig/03.AOD_distribution_validation/aeronet_vs_misr_eachyear_mean.png"),
         width = 7000,height = 2500,res=200)
print(distri_compare_each)
dev.off()

##---------** Distribution validation of each site----------
distri_valid_allyear <- readRDS(paste0("result/Rdata/Aerronet_AOD_valid/aerronet_vs_misr_daily_distribution.Rds"))%>%
  mutate(Site_adjust = ifelse(AERONET_Site%in%sites_sel,AERONET_Site,"Other sites"))%>%
  mutate(across(Site_adjust, factor, levels = site_order))%>%
  rename(MISR = VD)%>%
  rename(AERONET = VD_obs)

valid_ID <- levels(distri_valid_allyear$Site_adjust)
for (st in 1:length(valid_ID)) {
  distri_valid_st_allyear <- distri_valid_allyear%>%
    filter(Site_adjust == valid_ID[st])

  ##  distribution for each site
  distri_valid_st <- distri_valid_st_allyear%>%
    group_by(Lon,Lat,AERONET_Site,radius)%>%
    summarise(across(c("MISR","AERONET"), mean))%>%
    as.data.frame()%>%
    melt(measure.vars=c("MISR","AERONET"))
  
  ## logmathic size bin: equal distance
  label <- seq(0,15, by=15/21)
  distri_valid_st$label <- label
  
  st_distri_compare <- curve.plot.facet(distri_valid_st)+
    labs(x = "Radius",y = "Derived dV/dr")+
    theme(legend.key.width=unit(1, "cm"))+
    theme(legend.key.height=unit(1, "cm"))+##legend height
    theme(legend.spacing.y = unit(.8, 'cm'))+
    theme(legend.text=element_text(colour = "black", size = rel(0.8)))+ ##legend text size
    theme(legend.title = element_blank())+
    theme(plot.subtitle = element_text(colour = "black", size = rel(1.25),hjust = 0.5))+
    theme(axis.text=element_text(size=rel(1)),
          axis.title=element_text(size=rel(1)))
  
  CairoPNG(paste0("result/Fig/03.AOD_distribution_validation/aeronet_vs_misr_site_",valid_ID[st],".png"),
           width = 3000,height = 1400,res=300)
  print(st_distri_compare)
  dev.off()

##  distribution for each validated day
  for (year in 2004:2016) {
    distri_valid_st_allyear
    distri_valid_st_year <- distri_valid_st_allyear%>%
      filter(Year == year)
    st_doy <- unique(distri_valid_st_year$DOY)%>%
      sort()
    if(length(st_doy) > 0){
      for (i in 1:length(st_doy)) {
        distri_valid_st_doy <- distri_valid_st_year%>%
          filter(DOY == st_doy[i])%>%
          melt(measure.vars=c("MISR","AERONET"))
        label <- seq(0,15, by=15/21)
        distri_valid_st_doy$label <- label

        st_distri_compare <- curve.plot.facet(distri_valid_st_doy)+
          labs(x = "Radius",y = "Derived dV/dr")+
          theme(legend.key.width=unit(1, "cm"))+
          theme(legend.key.height=unit(1, "cm"))+##legend height
          theme(legend.spacing.y = unit(.8, 'cm'))+
          theme(legend.text=element_text(colour = "black", size = rel(0.8)))+ ##legend text size
          theme(legend.title=element_blank())+
          theme(plot.subtitle = element_text(colour = "black", size = rel(1.25),hjust = 0.5))+
          theme(axis.text=element_text(size=rel(1)),
                axis.title=element_text(size=rel(1)))
        if(dir.exists(paste0("result/Fig/AOD_distribution_validation/",valid_ID[st]))){
          CairoPNG(paste0("result/Fig/03.AOD_distribution_validation/",valid_ID[st],"/aeronet_vs_misr_site_",valid_ID[st],
                          "_Year_",year,"_doy ",st_doy[i],".png"),width = 3000,height = 1400,res=300)
          print(st_distri_compare)
          dev.off()
        }else{
          dir.create(paste0("result/Fig/03.AOD_distribution_validation/",valid_ID[st]))
          CairoPNG(paste0("result/Fig/03.AOD_distribution_validation/",valid_ID[st],"/aeronet_vs_misr_site_",valid_ID[st],
                          "_Year_",year,"_doy ",st_doy[i],".png"),width = 3000,height = 1400,res=300)
          print(st_distri_compare)
          dev.off()
        }
      }
    }##  end of doy
  }##  end of year
}
