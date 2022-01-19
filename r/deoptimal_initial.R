# use DEoptim to calculate a reasonable starting value for mcmc
# library(DEoptim)
library(foreach)
library(doParallel)
# Function to be passed to DE optim
model.de.func <- function(pars,dat,bucket.size,swc.in.wilt,swc.in.cap,day.lag,use.smooth,q.given,q.s.given){
  
  # deal with missing q
  if(is.null(q.given)){
    q.val = pars[5]
  }else{
    q.val = q.given
  }
  
  if(is.null(q.s.given)){
    q.s.val = pars[6]
  }else{
    q.s.val = q.given
  }
# make prediction
  hufken.pace.pred <- phenoGrass.func.v13(dat,
                                          f.h = 222,
                                          f.t.opt = pars[1],
                                          f.extract = pars[2],
                                          f.sec= pars[3],
                                          f.growth = pars[4],
                                          q =  q.val,
                                          q.s =  q.s.val,
                                          bucket.size = bucket.size,
                                          swc.wilt = swc.in.wilt ,
                                          swc.capacity = swc.in.cap ,
                                          t.max = 45,
                                          day.lay = day.lag,use.smooth = use.smooth)
  
  
  #standardise with sd
  sd.gcc <- sd(hufken.pace.pred$cover,na.rm = T)

  resid.gs <- ((hufken.pace.pred$cover.hufken - hufken.pace.pred$cover)/sd.gcc)^2
  
  return(sum(resid.gs,na.rm = T))
}



# func to use DEoptim to calaulate initial values
get.ini.func <- function(par.df,...){
  # on.exit(stopCluster(cl = NULL))
  # setting control parameters and limits to values
  lower <- unname(par.df['min',])
  upper <- unname(par.df['max',]) 
  NPmax <- 100
  maxiter <- 200
  # 
  set.seed(1935)
  OptBB.de.fit <- DEoptim(fn=model.de.func,lower=lower,upper=upper,
                          dat=gcc.met.pace.df.16,
                          DEoptim.control(VTR = 1,
                                          NP = NPmax,itermax=maxiter,trace=1,parallelType = 0,
                                          parVar = list('phenoGrass.func.v13',
                                                        'pet.func',
                                                        't.func',
                                                        'drainage.func','constants'),#globalenv(),
                                          packages=list("lubridate",'Evapotranspiration')),
                          use.smooth = use.smooth,
                          swc.in.cap = swc.capacity,
                          swc.in.wilt = swc.wilt,
                          bucket.size = bucket.size,
                          day.lag=day.lag,...)
  # Sys.sleep(10)
  initial.vec <- unname(OptBB.de.fit$optim$bestmem)
  return(initial.vec)
}


