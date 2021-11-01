library(Evapotranspiration)
# library(doBy)
# library(zoo)
# pet fucntion####
data("constants") #this is used for penman et value; data from Adleide
pet.func <- function(Date,PPFD,Tair,Tmax,Tmin,RHmax,RHmin,u2,P = 101.3,lat = 33.618891 ){
  gcc.met.df <- data.frame(Date = as.Date(Date),
                           PPFD = PPFD,
                           Tair = Tair,
                           Tmax = Tmax,
                           Tmin = Tmin,
                           RHmax = RHmax,
                           RHmin = RHmin,
                           u2=u2)
  
  
  gcc.met.df$J <- lubridate::yday(gcc.met.df$Date)
  gcc.met.df$n <-  gcc.met.df$PPFD * 10^-6/4.57 * 3600 *24
  R_s <-  gcc.met.df$PPFD * 10^-6/4.57 * 3600 *24
  # data <- gcc.met.df
  
  Ta <- gcc.met.df$Tair
  P <- 101.3 
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta + 237.3)))/((Ta + 237.3)^2)
  gamma <- 0.00163 * P/constants$lambda
  d_r2 <- 1 + 0.033 * cos(2 * pi/365 * gcc.met.df$J)
  delta2 <- 0.409 * sin(2 * pi/365 * gcc.met.df$J - 1.39)
  
  w_s <- acos(tan(lat * pi /180) * tan(delta2))
  N <- 24/pi * w_s
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * 
                                               sin(delta2) + cos(constants$lat_rad) * cos(delta2) * 
                                               sin(w_s))
  R_so <- (0.75 + (2 * 10^-5) * constants$Elev) * R_a
  
  # R_s <- (constants$as + constants$bs * (gcc.met.df$n/N)) * R_a
  
  vs_Tmax <- 0.6108 * exp(17.27 * gcc.met.df$Tmax/(gcc.met.df$Tmax + 
                                                     237.3))
  vs_Tmin <- 0.6108 * exp(17.27 * gcc.met.df$Tmin/(gcc.met.df$Tmin + 
                                                     237.3))
  vas <- (vs_Tmax + vs_Tmin)/2
  vabar <- (vs_Tmin * gcc.met.df$RHmax/100 + vs_Tmax * gcc.met.df$RHmin/100)/2
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * 
    ((gcc.met.df$Tmax + 273.2)^4 + (gcc.met.df$Tmin + 273.2)^4)/2 * 
    (1.35 * R_s/R_so - 0.35)
  
  alpha <- 0.67
  R_ns <- (1 - alpha) * R_s
  R_n = R_ns - R_nl
  # if (windfunction_ver == "1948") {
  f_u = 2.626 + 1.381 * gcc.met.df$u2
  #   # else if (windfunction_ver == "1956") {
  #     f_u = 1.313 + 1.381 * u2
  
  Ea = f_u * (vas - vabar)
  Epenman.Daily <- delta/(delta + gamma) * (R_n/constants$lambda) + 
    gamma/(delta + gamma) * Ea
  Epenman.Daily <- max(0,Epenman.Daily)
  return(Epenman.Daily)
}

# get a scaling factor####
scaling.f.func <- function(map,f.h){
  # map is input, h should be fitted
  map / (map + f.h)
}

# t response function####
t.func <- function(t.mean,f.t.opt,t.max){
  # return((t.max-t.mean)/(t.max-f.t.opt)*(t.mean/f.t.opt)^(f.t.opt/(t.max-f.t.opt)))

  h.val <- pnorm(t.mean,mean=f.t.opt,sd = 2,lower.tail = F)
  l.val <- pnorm(t.mean,mean=f.t.opt,sd = 2,lower.tail = T)

  return(min(h.val,l.val)*2)
}
# drainage
drainage.func <- function(theta,
                          theta.sat = 0.3,
                          sigma = 23.3,
                          k.sat = 0.15 * 1000#mm d-1
){
  if(theta.sat<0.3) theta.sat=0.3
  
  theta.frac <- theta/theta.sat
  theta.frac <- max(min(1,theta.frac),0)
  
  # taken from SGS, Johnson 2005 P71
  q = k.sat*(theta.frac)^sigma 
  
  # perscision control
  if(q<0.001){q=0}
  return(q)
}

