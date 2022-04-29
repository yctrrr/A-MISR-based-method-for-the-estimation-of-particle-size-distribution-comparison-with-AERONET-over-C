library(ggplot2)
##-----------Facet plot-----------
# data <- multiyear_mean_V_result
# label = "2011"
# value = "v_value"
# facet = "radius"
# x = "x"
# y = "y"
# facet_text =30

grid_facet_plot <- function(data = data, label = NULL, x = "x", y = "y", value = "value", facet = "factor",
                            facet_text = 30){
  ggplot()+
    geom_tile(data = data, aes_string(x=x, y=y, fill=value), alpha=1)+
    facet_wrap(facet, nrow = 1, labeller = labeller(groupwrap = label_wrap_gen(10)))+
    labs(x= label, y=NULL)+##label
    coord_equal()+
    #theme 
    theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black",size=1))+
    theme(strip.background = element_rect(
      color="black", fill="white", size = 1))+ ##change the facetwrap background
    theme(strip.text.x = element_text(
      size = facet_text, color = "black",margin = margin(.1, -.1, .1, -.1, "cm")))+
    theme(axis.title = element_text(size=rel(2)))+
    theme(legend.position="right") +
    theme(legend.key.width=unit(0.8, "cm"))+
    theme(legend.key.height=unit(2.5, "cm"))+##legend height
    theme(legend.text = element_text(colour = "black", size = rel(1.2)))+ ##legend text size
    theme(legend.title = element_text(colour = "black", size = rel(1.5),angle = -90))
}

##-----------Single plot-----------
# value = "discrete"
grid_plot <- function(data = data, label = NULL, x = "x", y = "y", value = "value"){
  ggplot()+
    geom_tile(data = data, aes_string(x = x, y = y, fill = value), alpha=1)+
    labs(x= label, y=NULL)+##label
    coord_equal()+
    theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black",size = 1))+##panel line width
    theme(strip.background = element_rect(color="black", fill="white", size = 1))+ ##change the facetwrap background
    theme(strip.text.x = element_text(size = 30, color = "black",margin = margin(.1, -.1, .1, -.1, "cm")))+
    theme(axis.title = element_text(size = rel(2)))+
    theme(legend.position = "right") +
    theme(legend.key.width = unit(1, "cm"))+
    theme(legend.key.height = unit(1, "cm"))+##legend height
    theme(legend.spacing.y = unit(.8, 'cm'))+
    theme(legend.text = element_text(colour = "black", size = rel(1.2)))+ ##legend text size
    theme(legend.title = element_text(colour = "black", size = rel(1.5),angle = -90))
}

##-------Add coordinates---------
scale_x_longitude <- function(xmin = -180, xmax = 180, step = 1, ...) {
  xbreaks <- seq(xmin,xmax,step)
  xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, paste0(x,"^o", "*W"), ifelse(x > 0, paste0(x,"^o", "*E"),x))))
  return(scale_x_continuous(breaks = xbreaks, labels = parse(text = xlabels), ...))
}

scale_y_latitude <- function(ymin = -90, ymax = 90, step = 0.5, ...) {
  ybreaks <- seq(ymin,ymax,step)
  ylabels <- unlist(lapply(ybreaks, function(x) ifelse(x < 0, paste0(x,"^o", "*S"), ifelse(x > 0, paste0(x,"^o", "*N"),x))))
  return(scale_y_continuous(breaks = ybreaks, labels = parse(text = ylabels), ...))
}

add_coordinates <- function(base_map, xmin = -180, xmax = 180, xstep = 1, ymin =-90, ymax = 90, ystep = 1, size = 2){
  base_map +
    scale_x_longitude(xmin = xmin, xmax = xmax, step = xstep)+ 
    scale_y_latitude(ymin = ymin, ymax = ymax, step = ystep) +
    theme(axis.text = element_text(size = rel(size),colour = "black", face = "bold"))
}

