######################################
# read and process the met, swc, and gcc for YM
######################################

# process ym data
source('r/pace_data_process.R')

get.paddock.func <- function(treat){
  # read swc
  flux.swc.df <- readRDS("cache/flux_twoer/flux_swc.df.rds")
  
  # get gcc for specific treat
  if(treat == 'control'){
    tmp.gcc.df <- readRDS("cache/flux_twoer/gcc_control.rds")
    tmp.gcc.df$SubplotID <- 1000
    p.con <- 'Control'
    tmp.gcc.df$Precipitation <- 'Control'
    # tmp.gcc.df$Date <- strptime(tmp.gcc.df$)
    
    swc.sub.df <- flux.swc.df[flux.swc.df$treat == 'control',]
  }
  if(treat == 'irrigated'){
    tmp.gcc.df <- readRDS("cache/flux_twoer/gcc_irrigated.rds")
    tmp.gcc.df$SubplotID <- 1000
    p.con <- 'Irrigated'
    tmp.gcc.df$Precipitation <- 'Irrigated'
    
    swc.sub.df <- flux.swc.df[flux.swc.df$treat == 'irrigated',]
  }
  #get daily swc 
  library(doBy)
  swc.sub.df$Date <- as.Date(swc.sub.df$TIMESTAMP)
  swc.sub.daily.df <- summaryBy(VWC_Avg~Date,data = swc.sub.df,
                                FUN=median,na.rm=T,keep.names = T)
  names(swc.sub.daily.df) <- c('Date','vwc')
  
  # 
  tmp.gcc.df$Temperature <- 'Ambient'
  tmp.gcc.df$Species <- 'flux'
  # 
  library(doBy)
  tmp.gcc.df.daily <- summaryBy(.~Date +SubplotID +Temperature +
                                  Precipitation+ Species,data = tmp.gcc.df,
                                FUN=mean,na.rm=T,keep.names = T)
  
  tmp.gcc.df.daily.sd <- summaryBy(GCC~Date +SubplotID +Temperature +
                                  Precipitation+ Species,data = tmp.gcc.df,
                                FUN=sd,na.rm=T,keep.names = F)
  
  tmp.gcc.df.daily$GCC.sd <- tmp.gcc.df.daily.sd$GCC.sd
  
  gcc.date <- range(tmp.gcc.df.daily$Date)
  
  if(p.con=='Irrigated'){
    gcc.date[1] <- as.Date('2020-1-7')
  }
  
  # get met
  ym.met.df <- readRDS('cache/ym/ym.met.rds')
  ym.met.df <- subset(ym.met.df,select = -c(sensor.no,position))
  ym.met.df <- ym.met.df[!duplicated(ym.met.df),]
  
  ym.met.df.daily <- summaryBy(.~Date,data = ym.met.df[,!(names(ym.met.df) == 'Rain_mm_Tot')],
                               FUN=mean,na.rm=T,keep.names = T)
  
  ym.met.df.daily.rain <- summaryBy(Rain_mm_Tot~Date,data = ym.met.df[ym.met.df$plot.no == ym.met.df$plot.no[1],],
                               FUN=sum,na.rm=T,keep.names = T)
  
  ym.met.df.daily$Rain <- ym.met.df.daily.rain$Rain_mm_Tot
  
  ym.met.flu.swc.df <- merge(ym.met.df.daily,swc.sub.daily.df,
                                              by='Date',all.x = T)
  
  # 
  gcc.met.ym.df <- merge(tmp.gcc.df.daily,ym.met.flu.swc.df,by=c('Date'),all=T)
  

  # make names to  be the same
  gcc.met.ym.df$Temperature <- 'Ambient'
  gcc.met.ym.df$Species <- 'flux'
  gcc.met.ym.df$SubplotID <- 1000
  gcc.met.ym.df$Precipitation <- p.con
  # gcc.met.ym.df$vwc <- gcc.met.ym.df$swc
  gcc.met.ym.df$PAR.ros <- gcc.met.ym.df$PAR
  gcc.met.ym.df$rh <- gcc.met.ym.df$RH
  gcc.met.ym.df$irrig.tot <- gcc.met.ym.df$Rain
  gcc.met.ym.df$WS_ms_Avg <- gcc.met.ym.df$u2
  
  gcc.met.ym.df$harvest = 0

  gcc.met.ym.df <- gcc.met.ym.df[gcc.met.ym.df$Date <gcc.date[2]& gcc.met.ym.df$Date > gcc.date[1],]
  
  return(gcc.met.ym.df)
}
