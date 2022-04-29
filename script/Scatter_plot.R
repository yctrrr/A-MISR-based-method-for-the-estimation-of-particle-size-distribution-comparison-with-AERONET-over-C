library(ggplot2)
library(egg)
library(dplyr)
library(tidyr)

# cv.df <- allyear_df
# cv.df <- distri_valid_allyear
# cex=3;alpha=0.2;bins=30;binwidth=c(0.05, 0.05)
# coord_ylim=0.95;ylimrange=0.05
# xlim=NULL;ylim=NULL
# nrow =2

cv.plot <- function(cv.df,count=FALSE,path=NULL,filename=NULL,
                    xlim=NULL,ylim=NULL,xlab=NULL,
                    ylab=NULL,cex=3,alpha=0.2,bins=30,binwidth=c(0.05, 0.05),
                    coord_ylim=0.95,ylimrange=0.05){
  coord_x <- min(xlim)
  coord_y <- max(ylim) - min(ylim)
  R <- cor(cv.df$y, cv.df$ypred)
  RMSE <- sqrt(mean((cv.df$y - cv.df$ypred)^2))
  MAPE <- mean(abs(cv.df$y - cv.df$ypred))
  N <- nrow(cv.df)
  model <- lm(cv.df$ypred ~ cv.df$y)
  c <- round(coef(model),3)
  if (c[1]<0){
    fm <- paste("Y=", c[2],"X",c[1],sep="")
  }else{
    fm <- paste("Y=", c[2],"X+",c[1],sep="")
  }
  ex=c(fm,      
       paste0("R=", formatC(round(R,digits = 4), digits=2,format="fg", flag="#")),
       paste0("MAE=",  formatC(round(MAPE,digits = 4), digits=3,format="fg", flag="#")),
       paste0("RMSE=", formatC(round(RMSE,digits = 4), digits=3,format="fg", flag="#")),
       paste0("N=",N)
  )

  cv.df %>%
    ggplot(aes(y,ypred)) +
    geom_point(size = 1.5, color = "#2b8cbe", alpha = 0.5) +
    annotate("text", x = coord_x*0.96 , y = c(coord_y*coord_ylim, coord_y*(coord_ylim-ylimrange),
                                            coord_y*(coord_ylim- ylimrange*2),coord_y*(coord_ylim-ylimrange*3),
                                            coord_y*(coord_ylim- ylimrange*4)),
             label = ex,hjust = 0, cex = cex)+
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme(panel.background = element_rect(size = 1,fill='white', colour='black'),
          panel.border = element_rect(size = 1,colour = "black", fill=NA))+
    geom_smooth(method = "lm",
                se = FALSE, colour="black",
                size = 0.5, fullrange=TRUE)+
    geom_abline(slope=1,linetype = "dashed",
                intercept=0, colour="blue", size=0.25)
  
}


cv.plot.facet.site <- function(cv.df,count=FALSE,xlim=NULL,ylim=NULL,path=NULL,filename=NULL,xlab=NULL,
                    ylab=NULL,nrow=NULL,cex=5,alpha=0.2,bins=30,binwidth=c(0.05, 0.05),
                    coord_ylim=0.95,ylimrange=0.05, vjust = 1.25){

  # site_name <- unique(cv.df$Site)
  site_name <- levels(cv.df$Site)
  for (i in 1:length(site_name)) {
    cv.df.sel <- filter(cv.df, Site == site_name[i])
    R <- cor(cv.df.sel$y, cv.df.sel$ypred)
    RMSE <- sqrt(mean((cv.df.sel$y - cv.df.sel$ypred)^2))
    MAPE <- mean(abs(cv.df.sel$y - cv.df.sel$ypred))
    N <- nrow(cv.df.sel)
    model <- lm(cv.df.sel$ypred ~ cv.df.sel$y)
    c <- round(coef(model),3)
    if (c[1]<0){
      fm <- paste("Y=", c[2],"X",c[1],sep="")
    }else{
      fm <- paste("Y=", c[2],"X+",c[1],sep="")
    }
    ex <- paste0(fm,"\nR=", formatC(round(R,digits = 4), digits=2, format="fg", flag="#"),
                 "\nMAE=", formatC(round(MAPE,digits = 4), digits=3, format="fg", flag="#"),
                 "\nRMSE=", formatC(round(RMSE,digits = 4), digits=3, format="fg", flag="#"),
                 "\nN=",N)
    
    if (i == 1){
      my_tag <- ex
    }else{
      my_tag <- c(my_tag, ex)
    }
  }
  cv.df


  facet_plot <- cv.df %>%
    ggplot(aes(y,ypred)) +
    geom_point(size = 1, aes(color = Site), alpha = 0.5,shape=20) +
    # facet_wrap(.~Site, nrow = nrow, scales = "free") +
    facet_wrap(.~Site, nrow = nrow) +
    coord_cartesian(xlim = xlim, ylim=ylim) +
    stat_smooth(method = "lm",
                se = FALSE,colour="black",
                size=0.1, fullrange = TRUE)+
    geom_abline(slope = 1,linetype = "dashed",
                intercept=0, colour="blue",size=0.1)+
    theme(panel.background = element_rect(size=0.5,fill='white', colour='black'),
          panel.border = element_rect(size=0.5,colour = "black", fill = NA))


    tag_facet(facet_plot,
              x = 0.01, y = Inf, 
              vjust = vjust, hjust = 0,
              open = "", close = "",
              fontface = 1,
              size = cex,
              tag_pool = my_tag)
    
}


## label:distinguish different facets
# cv.df <- distri_valid_allyear
# facet = "Year"
# label = "radius<"
# nrow = 1
# xlim = ylim = c(0,1)

cv.plot.facet.blue <- function(cv.df,count=FALSE,xlim=NULL,ylim=NULL,path=NULL,filename=NULL,xlab=NULL,
                               ylab=NULL, nrow=NULL, cex=5, alpha=0.2, bins=30, binwidth=c(0.05, 0.05),
                               facet = "Year", label = "Year=", coord_ylim=0.95,ylimrange=0.05, vjust = 1.25){
  
  facet_seq <- sort(unique(cv.df[,facet]))
  for (i in 1:length(facet_seq)) {
    cv.df.sel <- filter(cv.df, (!!as.name(facet) == facet_seq[i]))
    R <- cor(cv.df.sel$y, cv.df.sel$ypred)
    RMSE <- sqrt(mean((cv.df.sel$y - cv.df.sel$ypred)^2))
    MAPE <- mean(abs(cv.df.sel$y - cv.df.sel$ypred))
    N <- nrow(cv.df.sel)
    model <- lm(cv.df.sel$ypred ~ cv.df.sel$y)
    c <- round(coef(model),3)
    if (c[1]<0){
      fm <- paste("Y=", c[2],"X",c[1],sep="")
    }else{
      fm <- paste("Y=", c[2],"X+",c[1],sep="")
    }
    ex <- paste0(fm,"\nR=", formatC(round(R,digits = 4), digits=2, format="fg", flag="#"),
                 "\nMAE=", formatC(round(MAPE,digits = 4), digits=3, format="fg", flag="#"),
                 "\nRMSE=", formatC(round(RMSE,digits = 4), digits=3, format="fg", flag="#"),
                 "\nN=",N,
                 "\n",label,facet_seq[i])

    if (i == 1){
      my_tag <- ex
    }else{
      my_tag <- c(my_tag, ex)
    }
  }
  cv.df
  
  # facetlims <- cv.df%>%
  #   group_by(Year)%>%
  #   summarise(min = min(y, ypred), max = max(y, ypred))%>%
  #   gather(range, y, -Year)%>%
  #   mutate(ypred = y, range = NULL)%>%
  #   as.data.frame()

  facet_plot <- cv.df %>%
    ggplot(aes(y,ypred)) +
    geom_point(size = 1, color = "#2b8cbe", alpha = 0.5,shape=20) +
    # facet_wrap(.~Site, nrow = nrow, scales = "free") +
    facet_wrap(facet, nrow = nrow) +
    coord_cartesian(xlim = xlim, ylim=ylim) +
    stat_smooth(method = "lm",
                se = FALSE,colour="black",
                size=0.1, fullrange = TRUE)+
    geom_abline(slope = 1,linetype = "dashed",
                intercept=0, colour="blue",size=0.1)+
    theme(panel.background = element_rect(size=0.5,fill='white', colour='black'),
          panel.border = element_rect(size=0.5,colour = "black", fill = NA))
  
  tag_facet(facet_plot,
            x = 0.01, y = Inf, 
            vjust = vjust, hjust = 0,
            open = "", close = "",
            fontface = 1,
            size = cex,
            tag_pool = my_tag)
}


curve.plot.facet <- function(df,xlim=NULL,ylim=NULL,xlab=NULL,limits =NULL, labels = NULL,
                             ylab=NULL,cex=3,alpha=0.2,bins=30,binwidth=c(0.05, 0.05),
                             coord_ylim=0.95,ylimrange=0.05){
  ggplot(df,aes(label,value))+
    geom_point(aes(color = variable))+
    facet_wrap(.~variable,scales = "free")+ 
    stat_smooth(method = "gam",
                se = FALSE,colour="black",
                size=0.1, fullrange = TRUE)+
    theme(panel.background = element_rect(size=0.5,fill='white', colour='black'),
          panel.border = element_rect(size=0.5,colour = "black", fill = NA))
}
