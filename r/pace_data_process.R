# #######################################
# cleans and merges the pace data 
# met, gcc, swc
# #######################################

gcc.met.pace.df <- readRDS('cache/gcc.met.pace.df.rds')
gcc.met.pace.df = gcc.met.pace.df[gcc.met.pace.df$Date >= as.Date('2018-10-1')&
                                    gcc.met.pace.df$Date <= as.Date('2019-12-1'), ]
# with(gcc.met.pace.df[gcc.met.pace.df$SubplotID == 'S1P1A',],plot(GCC~Date))

# function to process GCC and met data
get.norm.gcc.func <- function(df,norm.min.max){
  if(is.null(norm.min.max)){
    quantiles.5.95 <- quantile(df$GCC[!is.na(df$GCC)],
                               c(.01,.99),na.rm=T)

    quantiles.5.95[1] <- 0.3
    # quantiles.5.95 <- c(0.3,0.43)
  }else{
    quantiles.5.95 <- norm.min.max
  }
  
  # quantiles.5.95[1] = 0.3197
  df$GCC.norm <- (df$GCC - quantiles.5.95[1]) /
    (quantiles.5.95[2] - quantiles.5.95[1])
  
  return(df)
}

# get smooth GCC with GAM
get.smooth.gcc.func = function(Date.vec,gcc.vec){
  library(mgcv)
  library(lubridate)
 
  gam.out.df = data.frame(x = seq_along((Date.vec)),
                          y = gcc.vec)
 
  gam.in.df <- gam.out.df#gam.out.df[!is.na(gam.out.df$y),]
  
  gam.frdm = round(length(gam.in.df$x)/3)
  
  fit.gam <- gam(y~s(x,k = gam.frdm),data = gam.in.df)
  
  gam.in.df$fitted.val = predict(fit.gam,gam.in.df)

  return(gam.in.df$fitted.val)
}

# get pace data by species or treatment
get.pace.func <- function(gcc.met.pace.df,
                          species.in,
                          prep.in,
                          temp.in,subplot = NA,
                          norm.min.max=NULL){
  # gather data by plot or treatment
  if(is.na(subplot)){
    temp.df <- gcc.met.pace.df[gcc.met.pace.df$Species == species.in&
                                 gcc.met.pace.df$Precipitation == prep.in&
                                 gcc.met.pace.df$Temperature == temp.in,]
  }else{
    temp.df <- gcc.met.pace.df[gcc.met.pace.df$SubplotID == subplot,]
  }

  temp.df <- temp.df[!is.na(temp.df$Date),]
  
  if(nrow(temp.df)<2)stop('species/treament incorrect')

  # get normalised GCC by plot
  subplot.vec <- unique(temp.df$SubplotID)
  
  temp.ls <- split(temp.df,temp.df$SubplotID)
  
  result.ls <- lapply(temp.ls,get.norm.gcc.func,norm.min.max=norm.min.max)
  
  temp.norm.df <- do.call(rbind,result.ls)
  # get sd and mean
  sd.df <- summaryBy(GCC.norm + vwc ~ Date,
                     data = temp.norm.df,FUN=c(sd),na.rm=TRUE,keep.names = F)
  
  mean.df <- summaryBy(.~ Date,
                       data = temp.norm.df,FUN=c(mean),na.rm=TRUE,keep.names = T)

  test.df <- merge(mean.df,sd.df)
  # remove missing values 
  start.date <- test.df$Date[min(which(!is.na(test.df$GCC.norm)))]
  
  gcc.met.pace.df.na.rm <- test.df[test.df$Date > start.date,]

  gcc.met.pace.df.16 <- rbind(gcc.met.pace.df.na.rm,
                              test.df[test.df$Date <= start.date &
                                        test.df$Date > (start.date- (day.lag + 1)),])
  
  gcc.met.pace.df.16 <- gcc.met.pace.df.16[order(gcc.met.pace.df.16$Date),]
  
  gcc.met.pace.df.16 <- gcc.met.pace.df.16[!is.na(gcc.met.pace.df.16$Date),]

   # 
  gcc.met.pace.df.16$irrig.tot[is.na(gcc.met.pace.df.16$irrig.tot)] <- 0

  # 
  gcc.met.pace.df.16$PPFD <- na.locf(gcc.met.pace.df.16$PAR.ros)
  gcc.met.pace.df.16$u2 <- na.locf(gcc.met.pace.df.16$WS_ms_Avg)
  gcc.met.pace.df.16$Rain <- na.locf(gcc.met.pace.df.16$irrig.tot)
  gcc.met.pace.df.16$Tair <- na.locf(gcc.met.pace.df.16$Tair)
  gcc.met.pace.df.16$RHmax <- na.locf(gcc.met.pace.df.16$RHmax)
  gcc.met.pace.df.16$RHmin <- na.locf(gcc.met.pace.df.16$RHmin)
  gcc.met.pace.df.16$Tmax <- na.locf(gcc.met.pace.df.16$Tmax)
  gcc.met.pace.df.16$Tmin <- na.locf(gcc.met.pace.df.16$Tmin)
  gcc.met.pace.df.16$rh <- na.locf(gcc.met.pace.df.16$rh)
  
  # set plot with no vwv to be 0.08
  # this allows the model to a a vwc to start with 
  #  doesn't do anything else
  if((sum(gcc.met.pace.df.16$vwc,na.rm=T))==0){
    gcc.met.pace.df.16$vwc = 0.08
    warning('VWC not given')
  }
  
  # use 10% of value as sd if sd not given
  # only used when 1 obs per day
  if(sum(gcc.met.pace.df.16$GCC.norm.sd,na.rm=T) == 0){
    # gcc.met.pace.df.16$GCC.norm.sd = abs( gcc.met.pace.df.16$GCC.norm * 0.1)
    gcc.met.pace.df.16$GCC.norm.sd = gcc.met.pace.df.16$GCC.norm.sd /(0.42-0.3)
    warning('SD of GCC not given')
  }

  gcc.met.pace.df.16 <- gcc.met.pace.df.16[order(gcc.met.pace.df.16$Date),]
  gcc.met.pace.df.16$GCC.norm.smooth = get.smooth.gcc.func(Date.vec = gcc.met.pace.df.16$Date, 
                                                           gcc.vec = gcc.met.pace.df.16$GCC.norm)
  

  gcc.met.pace.df.16$GCC.norm.smooth[gcc.met.pace.df.16$GCC.norm.smooth<0] <- 0
  gcc.met.pace.df.16$GCC.norm.smooth[gcc.met.pace.df.16$GCC.norm.smooth>1] <- 1

  
  if(is.null(gcc.met.pace.df.16$harvest)){
    gcc.met.pace.df.16$harvest = 0
  }
  
  print(paste0(species.in,' data processed'))
  return(gcc.met.pace.df.16)
}

