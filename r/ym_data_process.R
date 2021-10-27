# process ym data
source('r/pace_data_process.R')


get.ym.func <- function(treat){
  
  # c.wd <- getwd()
  # setwd('c:/repo/dn_gcc/')
  # on.exit(setwd(c.wd))
   # get gcc for specific treat
  if(treat == 'Control'){
    df.45 <- readRDS('cache/ym/gcc_45.rds')
    df.45$SubplotID <- 45
    df.18 <- readRDS('cache/ym/gcc_18.rds')
    df.18$SubplotID <- 18
    
    # 
    f.rain.reduce <- 1
    
    tmp.gcc.df <- rbind(df.18,df.45)
    tmp.gcc.df$Precipitation <- 'Control'
  }
  if(treat == 'Drought'){
    df.14 <- readRDS('cache/ym/gcc_14.rds')
    df.14$SubplotID <- 14
    df.27 <- readRDS('cache/ym/gcc_27.rds')
    df.27$SubplotID <- 27
    df.38 <- readRDS('cache/ym/gcc_38.rds')
    df.38$SubplotID <- 38
    
    f.rain.reduce <- 0.35
    
    tmp.gcc.df <-do.call(rbind,list(df.14,df.27,df.38))
    tmp.gcc.df$Precipitation <- 'Drought'
  }
  if(treat == 'up'){
    tmp.gcc.df <- readRDS('cache/ym/gcc_up.rds')
    tmp.gcc.df$Precipitation <- 'up'
    tmp.gcc.df$SubplotID <- 'up'
    
    f.rain.reduce <- 1
  }
  if(is.numeric(treat)){
    fn <- sprintf('cache/ym/gcc_%s.rds',treat)
    tmp.gcc.df <- readRDS(fn)
    if(treat %in% c(14,27,38)){
      prep <- 'Drought'
      f.rain.reduce <- 0.35
    }else{
      prep <- 'Control'
      f.rain.reduce <- 1
    }
    tmp.gcc.df$Precipitation <- prep
    
    tmp.gcc.df$SubplotID <- treat
  }
  

  tmp.gcc.df$Temperature <- 'Ambient'
  tmp.gcc.df$Species <- 'ym'
  
  library(doBy)
  tmp.gcc.df.daily <- summaryBy(.~Date +SubplotID +Temperature +
                            Precipitation+ Species,data = tmp.gcc.df,
                            FUN=mean,na.rm=T,keep.names = T)
  
  # get met
  ym.met.df <- readRDS('cache/ym/ym.met.rds')
  ym.met.df <- subset(ym.met.df,select = -c(sensor.no,position))
  ym.met.df <- ym.met.df[!duplicated(ym.met.df),]
  
  ym.met.df$SubplotID <- ym.met.df$plot.no
  gcc.met.ym.df <- merge(tmp.gcc.df.daily,ym.met.df,by=c('Date','SubplotID'))
  
  # get names to  be the same
  gcc.met.ym.df$vwc <- gcc.met.ym.df$swc
  gcc.met.ym.df$PAR.ros <- gcc.met.ym.df$PAR
  gcc.met.ym.df$rh <- gcc.met.ym.df$RH
  gcc.met.ym.df$irrig.tot <- gcc.met.ym.df$Rain * f.rain.reduce
  gcc.met.ym.df$WS_ms_Avg <- gcc.met.ym.df$u2
  
  gcc.met.ym.df$harvest = 0
  # # use the pace function to get standard gcc and met
  # gcc.met.ym.processed.df <- get.pace.func(gcc.met.ym.df,
  #                                  species.in = 'ym',
  #                                  prep.in = treat,
  #                                  temp.in = 'Ambient',
  #                                  subplot = NA)
  
  return(gcc.met.ym.df)
}
