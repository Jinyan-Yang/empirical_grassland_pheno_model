
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

  # # check whether it rained in the past days
  # rained.vec <- c()
  # rained.vec[1:day.lay] <- 0
  # for(i in (day.lay + 1):nrow(gcc.df)){
  #   if(sum(gcc.df$Rain[(i-day.lay):i],na.rm=T)>0){
  #     rained.vec[i] <- sum(gcc.df$Rain[(i-day.lay):i],na.rm=T)
  #   }else{
  #     rained.vec[i] <- 0
  #   }
  # }

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
   
    # default no growth & no sen
    g = 1
    d = 1

    # # grwoth after rain
    # if(rained.vec[nm.day] > 0){
    #   g = 1
    #   d = 0
    # }else{
    #   g = 1
    #   d = 1
    # }

    # # # growth if soil is really wet
    # if(!is.na(loss.f)){
    #   if(loss.f > 0.7){
    #     g = 1
    #     d = 0
    #   }
    # }

      # if(swc.norm > 0.7){
      #   g = 1
      #   d = 0
      # }
      # 
    
    # dormaince 
    # dor = gcc.df$rad.norm[nm.day]
    #min(max(dor,1),0)
    dor = 1
    
    # temperature effect
    g.value <-  t.func(t.m[nm.day],f.t.opt,t.max)

    # plant cover
    growth.vec[nm.day] <- dor*g * g.value * f.growth *
      loss.f *
      #water.avi.norm*
      (1 - cover.pred.vec[nm.day-1] / cover.max)
    
    senescence.vec[nm.day] <- d * f.sec * 
      loss.f.s *
      # 4*(1 - cover.pred.vec[nm.day-1])*
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
      # loss.f^2*
      loss.f.soil*
      # ((swc.vec[nm.day-1]/bucket.size - swc.wilt)/(swc.capacity-swc.wilt))^2 *
      et[nm.day]
    transp.vec[nm.day] <- f.extract *#g.value*
      swc.vec[nm.day-1] *
      # loss.f.soil*
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


# phenoGrass.func.v13 <- function(gcc.df,
#                                 f.h,
#                                 f.t.opt,
#                                 f.extract,
#                                 f.sec,
#                                 f.growth,
#                                 swc.wilt ,
#                                 swc.capacity ,
#                                 bucket.size,
#                                 t.max,
#                                 day.lay,
#                                 use.smooth=FALSE,
#                                 q = 1){
#   
#   if(f.t.opt < 0){
#     f.t.opt = 0
#   }
#   if(f.extract < 0){
#     f.extract = 0
#   }
#   if(f.sec < 0){
#     f.sec = 0
#   }
#   if(f.growth < 0){
#     f.growth = 0
#   }
#   if(q < 0){
#     q = 0
#   }
#   
#   # set the lag factor; in num of days
#   start.date <- gcc.df$Date[min(which(!is.na(gcc.df$GCC.norm)))]
#   
#   gcc.df <- gcc.df[gcc.df$Date > (start.date - day.lay),]
#   sf.value <- scaling.f.func(mean(gcc.df$map,na.rm=TRUE),f.h)
#   
#   # decide whether to use smooth gcc
#   if(use.smooth==TRUE){
#     gcc.df$cover <-  gcc.df$GCC.norm.smooth * sf.value
#   }else{
#     gcc.df$cover <-  gcc.df$GCC.norm * sf.value
#   }
#   
#   # set up the inital conditions
#   swc.vec <- gcc.df$vwc * bucket.size
#   et <- c()
#   cover.pred.vec <- c()
#   cover.pred.vec[day.lay] <- gcc.df$cover[!is.na(gcc.df$cover)][1]
#   water.avi <- c()
#   water.avi <- swc.vec 
#   water.lag <- c()
#   water.lag <- water.avi
#   t.m <- growth.vec <- senescence.vec <- evap.vec <- transp.vec <- c()
#   
#   # # check whether it rained in the past days
#   # rained.vec <- c()
#   # rained.vec[1:day.lay] <- 0
#   # for(i in (day.lay + 1):nrow(gcc.df)){
#   #   if(sum(gcc.df$Rain[(i-day.lay):i],na.rm=T)>0){
#   #     rained.vec[i] <- sum(gcc.df$Rain[(i-day.lay):i],na.rm=T)
#   #   }else{
#   #     rained.vec[i] <- 0
#   #   }
#   # }
#   
#   # calcualte the par values
#   cover.max <- max(gcc.df$cover,na.rm=TRUE)
#   # rad.min <-  min(gcc.df$PPFD,na.rm=TRUE)
#   # rad.max <-  max(gcc.df$PPFD,na.rm=TRUE)
#   # gcc.df$rad.norm <- (gcc.df$PPFD - rad.min) / (rad.max - rad.min)
#   
#   # model start
#   for (nm.day in (day.lay+1):nrow(gcc.df)){
#     
#     water.avi[nm.day] <- max(0,(swc.vec[nm.day-1]- swc.wilt*bucket.size))
#     
#     # define water stress using a beta function
#     swc.norm <- (water.avi[nm.day] / (swc.capacity - swc.wilt)/bucket.size)
#     loss.f <- swc.norm^q 
#     loss.f <- min(1,loss.f)
#     # assuming soil stress is the same 
#     loss.f.soil<-loss.f
#     
#     # # define the legency effect 
#     i=0
#     while(i+1<day.lay & (nm.day-i)>0){
#       i=i+1
#     }
#     
#     t.m[nm.day] <- gcc.df$Tair[nm.day]#mean(gcc.df$Tair[(nm.day-15):nm.day],na.rm=TRUE) #hufkens used 15 days
#     # hufkens used evaportanspiration from Hargreaves 1985
#     # here is from evapotranspiration R package
#     et[nm.day] <- pet.func(gcc.df$Date[nm.day],gcc.df$PPFD[nm.day],
#                            gcc.df$Tair[nm.day],gcc.df$Tmax[nm.day], gcc.df$Tmin[nm.day],
#                            gcc.df$RHmax[nm.day],gcc.df$RHmin[nm.day], gcc.df$u2[nm.day])    
#     g=1
#     if(water.lag[nm.day] <=  water.lag[nm.day-1]){
#       d = 1
#     }else{
#       d = 0
#     }
#     
#     # plant cover
#     g.value <- 1# t.func(t.m[nm.day],f.t.opt,t.max)
#     growth.vec[nm.day] <- g * g.value * f.growth * 
#       loss.f * 
#       #water.avi.norm*
#       (1 - cover.pred.vec[nm.day-1] / cover.max)
#     senescence.vec[nm.day] <- d *# f.sec * (1 - loss.f) * 
#       (1 - cover.pred.vec[nm.day-1] )*cover.pred.vec[nm.day-1]
#     
#     cover.pred.vec[nm.day] <- cover.pred.vec[nm.day-1] + growth.vec[nm.day] - senescence.vec[nm.day]                 
#     
#     # give min
#     cover.pred.vec[nm.day] <- max(0,min(cover.pred.vec[nm.day],cover.max))
#     
#     # account for harvest
#     if(gcc.df$harvest[nm.day] == 1){
#       cover.pred.vec[nm.day] <- gcc.df$cover[nm.day+2]
#     }
#     # calculate swc
#     evap.vec[nm.day] <- (1 - cover.pred.vec[nm.day-1]) * 
#       # loss.f^2*
#       loss.f.soil*#^2*
#       # ((swc.vec[nm.day-1]/bucket.size - swc.wilt)/(swc.capacity-swc.wilt))^2 *
#       et[nm.day]
#     transp.vec[nm.day] <- f.extract * 
#       # swc.vec[nm.day-1] * 
#       # water.avi.norm*
#       loss.f *
#       cover.pred.vec[nm.day] / cover.max
#     
#     swc.vec[nm.day] <- swc.vec[nm.day-1] + gcc.df$Rain[nm.day] - evap.vec[nm.day] - transp.vec[nm.day]
#     
#     swc.vec[nm.day] <- max(0,min(swc.capacity * bucket.size,swc.vec[nm.day]))
#     
#   }
#   gcc.df$ppt <- gcc.df$Rain
#   gcc.df$cover.hufken <- cover.pred.vec
#   gcc.df$swc.hufken <- swc.vec
#   gcc.df$vwc.hufken <- swc.vec / bucket.size
#   gcc.df$water.avi <- water.avi
#   gcc.df$evap <- evap.vec 
#   gcc.df$tran <- transp.vec
#   # out.df <- data.frame(gcc.df)
#   # out.df <- out.df[!is.na(out.df$cover),]
#   # print('model worked')
#   return(gcc.df)
# }
# 
