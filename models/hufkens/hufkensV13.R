
####################################################################################################################
##############################################V13 of the hufkens model##############################################
####################################################################################################################

# load need func and packages
source('models/hufkens/hufkens_common_fun.R')

# model function
phenoGrass.func.v13 <- function(gcc.df,
                                f.h,
                                f.t.opt,
                                f.extract,
                                f.sec,
                                f.growth,
                                swc.wilt ,
                                swc.capacity ,
                                bucket.size,
                                t.max,
                                day.lay,
                                use.smooth=FALSE,
                                q = 1,q.s=1){
  
  ####################################################################################################################
  # inputs
         # gcc.df: data frame with GCC, and met inputs
         # f.h: gcc to cover scaling factor; mm/yr
         # f.t.opt: optimal temperature for growth; degree C
         # f.extract: max transpiration rate; mm d-1 
         # f.sec & f.growth: intrinsic senescence and growth rate; cover d-1 
         # swc.wilt & swc.capacity: soil wilting point and capacity; decimal
         # bucket.size: rooting depth; mm
         # t.max; max T after which no growth allowed; degree C
         # day.lay; number of days plant wait to reponsd to rainfall
         # use.smooth; switch to use GAM smoothed GCC; TRUE- use smoothed data
         # q & q.s: power of the beta function for growth and senescence; unitless
  # outputs
        # data frame with met and predicted cover, SWC, runoff, drainage, evaporation
  ####################################################################################################################

  # model begine####
  # 1.set up initial conditions
  # ignore the dates without GCC
  start.date <- gcc.df$Date[min(which(!is.na(gcc.df$GCC.norm)))]

  gcc.df <- gcc.df[gcc.df$Date > (start.date - day.lay),]
  
  # MAP scaling factor
  sf.value <- 1#scaling.f.func(mean(gcc.df$map,na.rm=TRUE),f.h)

  # decide whether to use smooth gcc
  if(use.smooth==TRUE){
    gcc.df$cover <-  gcc.df$GCC.norm.smooth * sf.value
  }else{
    gcc.df$cover <-  gcc.df$GCC.norm * sf.value
  }
#  use drainage function to determine the real capacity
  loss.at.full <- drainage.func(swc.capacity,theta.sat= swc.capacity)
  swc.cap.true <- swc.capacity  #- round(loss.at.full,digits=2)/ bucket.size
  # set up the inital conditions
  if(length(unique(gcc.df$vwc)) == 1){
    # cover.fraction <- gcc.df$cover[!is.na(gcc.df$cover)][1]
    swc.vec <- (swc.cap.true - swc.wilt)* bucket.size * gcc.df$cover + swc.wilt*bucket.size
    # swc.vec[swc.vec<swc.wilt*bucket.size] <- swc.wilt*bucket.size
    # print(cover.fraction)
  }else{
    swc.vec <- gcc.df$vwc * bucket.size
  }
  
  swc.vec[day.lay] <- swc.vec[!is.na(swc.vec)][1]
 
  et <- c()
  cover.pred.vec <- c()
  cover.pred.vec[day.lay] <- gcc.df$cover[!is.na(gcc.df$cover)][1]
  water.avi <- c()
  water.avi <- swc.vec - swc.wilt*bucket.size
  water.lag <- c()
  water.lag <- water.avi
  t.m <- growth.vec <- senescence.vec <- evap.vec <- transp.vec <- runOff <-drain.vec <-  c()

  # calcualte the par values
  cover.max <- 1#max(gcc.df$cover,na.rm=TRUE)
  # rad.min <-  min(gcc.df$PPFD,na.rm=TRUE)
  # rad.max <-  max(gcc.df$PPFD,na.rm=TRUE)
  gcc.df$rad.norm <- 1#(gcc.df$PPFD - rad.min) / (rad.max - rad.min)

  # 2.daily simulation start
  for (nm.day in (day.lay+1):nrow(gcc.df)){
    water.avi[nm.day] <- max(0,(swc.vec[nm.day-1]- swc.wilt*bucket.size))

    # define water stress using a beta function
    swc.norm <-water.avi[nm.day] / (swc.cap.true - swc.wilt) / bucket.size
    swc.norm <- max(0,min(1,swc.norm))
    loss.f <- swc.norm^q
    loss.f <- min(1,loss.f)
    # assuming sene stress is the not same
    loss.f.s <- (1-swc.norm)^q.s
    loss.f.s <- min(1,loss.f.s)
    # assume soil evap is linear to swc
    loss.f.soil <- swc.norm
    
    # # # define the legency effect
    # i=0
    # while(i+1<day.lay & (nm.day-i)>0){
    #   i=i+1
    # }

    days.past <- max(c(1,(nm.day-15))) #hufkens used 15 days
    t.m[nm.day] <- gcc.df$Tmax[nm.day]#mean(gcc.df$Tmax[days.past:nm.day],na.rm=TRUE) #
    # hufkens used evaportanspiration from Hargreaves 1985
    # here is from evapotranspiration R package
    et[nm.day] <- pet.func(gcc.df$Date[nm.day],gcc.df$PPFD[nm.day],
                           gcc.df$Tair[nm.day],gcc.df$Tmax[nm.day], gcc.df$Tmin[nm.day],
                           gcc.df$RHmax[nm.day],gcc.df$RHmin[nm.day], gcc.df$u2[nm.day])    
   
    # default  growth &  sen#these control parameters
    # not used for now but kept for future use
    g = 1
    d = 1
    dor = 1
    
    # temperature effect
    g.value <-  t.func(t.m[nm.day],f.t.opt,t.max)

    # plant cover
    growth.vec[nm.day] <- dor*g * g.value * f.growth *
      loss.f *
      (1 - cover.pred.vec[nm.day-1] / cover.max)
    
    senescence.vec[nm.day] <- d * f.sec * 
      loss.f.s *
      cover.pred.vec[nm.day-1]/ cover.max

    cover.pred.vec[nm.day] <- cover.pred.vec[nm.day-1] + growth.vec[nm.day] - senescence.vec[nm.day]

    # give min
    cover.pred.vec[nm.day] <- max(0,min(cover.pred.vec[nm.day],cover.max))

    # account for harvest
    if(gcc.df$harvest[nm.day] == 1){
      cover.pred.vec[nm.day] <- mean(gcc.df$cover[nm.day:(nm.day+2)])
    }
    # calculate swc
    evap.vec[nm.day] <- (1 - cover.pred.vec[nm.day-1]) *
      loss.f.soil*
      et[nm.day]
    
    transp.vec[nm.day] <- f.extract *#g.value* Not used T dependence on ET

      loss.f.soil*
      cover.pred.vec[nm.day] / cover.max

    swc.vec[nm.day] <- swc.vec[nm.day-1] + gcc.df$Rain[nm.day] - evap.vec[nm.day] - transp.vec[nm.day]

    # assuming runoff  = whatever is over the capacity
    runOff[nm.day] <- max(0,(swc.vec[nm.day] - swc.capacity * bucket.size))
    
    swc.vec[nm.day] <- min(swc.capacity * bucket.size,swc.vec[nm.day])
    
    # apply drainage
    drain.vec[nm.day] <- drainage.func(swc.vec[nm.day] / bucket.size,
                                       theta.sat  = swc.capacity)
    swc.vec[nm.day] <- swc.vec[nm.day] - drain.vec[nm.day]

    swc.vec[nm.day] <- max(0,swc.vec[nm.day])
  }
  
  # 3.organise output
  gcc.df$ppt <- gcc.df$Rain
  gcc.df$cover.hufken <- cover.pred.vec
  gcc.df$swc.hufken <- swc.vec
  gcc.df$vwc.hufken <- swc.vec / bucket.size
  gcc.df$water.avi <- water.avi
  gcc.df$evap <- evap.vec
  gcc.df$tran <- transp.vec
  gcc.df$pet <- et
  gcc.df$runOff <- runOff
  gcc.df$drainage <- drain.vec
  gcc.df$growth <- growth.vec
  gcc.df$senescence <- senescence.vec
  
  return(gcc.df)
}

