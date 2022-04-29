##  load MISR AOD ratio
library(zoo)
##  function to reconstruct data frame
framer <- function(data){
  return (data.frame(y=data$y, x=data$x))
}

integral <- function(density,lower,upper,n){
  sel <- density$x <= upper&density$x >= lower
  integral <- sum(rollmean(density$x[sel],2)^n*diff(density$x[sel])*rollmean(density$y[sel],2))
  return(integral)
}

#------V3: calculate size distribution and volume of particle with given radius(aeronet radius)-------
## aod: MISR AOD with components values
## rm_valid: radius of interest
misr_aod_to_distri <- function(aod, rm_valid = c(0.050000,0.065604,0.086077,0.112939,0.148184,0.194429,0.255105,0.334716,0.439173,0.5,0.576227,0.756052,
                                                    0.991996,1.25,1.301571,1.707757,2,2.240702,2.939966,3.857452,5,5.061260,6.640745,8.713145,11.432287,15.000000,50)){
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

  ##  Number-weighted mode radius of 8 components
  radius_list <- c(0.03,0.06,0.12,1,0.06,0.06,0.5,1)

  number <- c(1,2,3,6,8,14,19,21)

  ext_list <- c(0.00039594,0.013401464,0.18231007,16.189339,
                0.014069699,0.014874,3.168676,15.51)
  
  ##  Component density calculation
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

  # integral(d1,0,100,0)
  for (comp_id in 1:8) {
    component <- aod%>%
      as.data.frame()%>%
      .[,paste0("Component_",formatC(number[comp_id],width = 2,flag=0))]
    
    ext <- ext_list[comp_id]
    range <- get(paste0("range",number[comp_id]))
    d <- get(paste0("d",number[comp_id]))
    ##  calculation
    integral_square <- integral(d,range[1],range[2],2)
    N <- component/(integral_square*3.14)
    
    for (r in 1:length(rm_valid)) {
      r_sel <- which.min(abs(d$x - rm_valid[r]))
      if(r==1){
        ## calculate normalized dN/dr and dV/dr
        N_distri_valid <- d[r_sel,]
        V_distri_valid <- mutate(N_distri_valid,y = y*3.14*4*x^3/3)
        
        ## calculate normalized particle 'number' and 'volume' for particle radius less than rm_valid um
        N_valid <- data.frame(x = rm_valid[r],
                              number = integral(d,0,rm_valid[r],0), component = d[r_sel, "component"])
        
        V_valid <- data.frame(x = rm_valid[r],
                              volume = 3.14*4*integral(d,0,rm_valid[r],3)/3, component = d[r_sel, "component"])
      }else{
        N_distri_valid <- full_join(d[r_sel,], N_distri_valid)
        V_distri_valid <- full_join(mutate(N_distri_valid,y = y*3.14*4*x^3/3), V_distri_valid)
        N_valid <- full_join(data.frame(data.frame(x = rm_valid[r],
                                                   number = integral(d,0,rm_valid[r],0), component = d[r_sel, "component"])),
                             N_valid)
        V_valid <- full_join(data.frame(x = rm_valid[r],
                                        volume = 3.14*4*integral(d,0,rm_valid[r],3)/3, component = d[r_sel, "component"]),
                             V_valid)
      }
    }##  end of valid radius

    ##  order the radius value
    N_distri_valid <- N_distri_valid[order(N_distri_valid$x),]
    V_distri_valid <- V_distri_valid[order(V_distri_valid$x),]
    N_valid <- N_valid[order(N_valid$x),]
    V_valid <- V_valid[order(V_valid$x),]
    
    ##  calculate dN/dr and dV/dr
    N_distri_valid_comp <- lapply(1:nrow(N_distri_valid), function(i){N*N_distri_valid$y[i]})
    V_distri_valid_comp <- lapply(1:nrow(V_distri_valid), function(i){N*V_distri_valid$y[i]})
    ##  calculate number and volume
    N_valid_comp <- lapply(1:nrow(N_valid),function(i){N*N_valid$number[i]})
    V_valid_comp <- lapply(1:nrow(V_valid),function(i){N*V_valid$volume[i]})
    
    ##  calculate the sum of number and volume for all components
    if(comp_id == 1){
      N_distri_join <- N_distri_valid_comp
      V_distri_join <- V_distri_valid_comp
      N_join <- N_valid_comp
      V_join <- V_valid_comp
    }else{
      N_distri_join <- lapply(seq_along(N_distri_valid_comp),function(i) N_distri_valid_comp[[i]]+N_distri_join[[i]])
      V_distri_join <- lapply(seq_along(V_distri_valid_comp),function(i) V_distri_valid_comp[[i]]+V_distri_join[[i]])
      N_join <- lapply(seq_along(N_valid_comp),function(i) N_valid_comp[[i]]+N_join[[i]])
      V_join <- lapply(seq_along(V_valid_comp),function(i) V_valid_comp[[i]]+V_join[[i]])
    }
  }## end of components


  # for (r in 1:length(rm_valid)) {
  #   ## match size distribution with AOD data
  #   aod <- mutate(aod,!!paste0("ND",rm_valid[r]):= N_distri_join[[r]])
  #   aod <- mutate(aod,!!paste0("VD",rm_valid[r]):= V_distri_join[[r]])
  #   aod <- mutate(aod,!!paste0("N",rm_valid[r]):= N_join[[r]])
  #   aod <- mutate(aod,!!paste0("V",rm_valid[r]):= V_join[[r]])
  # }
  # return(aod)

  for (r in 1:length(rm_valid)) {
    aod_rm <- mutate(aod, radius = rm_valid[r], ND = N_distri_join[[r]], VD = V_distri_join[[r]],
                     N = N_join[[r]], V = V_join[[r]])
    if(r == 1){
      aod_join <- aod_rm
    }else{
      aod_join <- full_join(aod_join, aod_rm)
    }
  }
  return(aod_join)
}
