# Fit model to all data sets#################################
# this file reads data from PACE, DN, and Flux tower
# and fit CH model to them
#############################################################

#read in data #################
ym.18.df <- get.ym.func(18)
gcc.met.con.df <- get.paddock.func('control')

# loop through all spcies/site##################
for (i in seq_along(species.vec)){
  
  # use different soil water cap and wilt for different site
  if(species.vec[i]=='ym'){
    df = ym.18.df
    ym.met.df <- readRDS('cache/ym/ym.met.rds')
    
    swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    swc.cap = round(swc.ym.con[[2]]*10)/10
    swc.wilt = round(swc.ym.con[[1]]*100)/100
    bucket.size=1000
  }else if(species.vec[i]=='flux'){
    df = gcc.met.con.df
    swc.q.con <- quantile(gcc.met.con.df$vwc,na.rm=T,probs = c(0.01,0.99))
    swc.cap =  round(swc.q.con[[2]]*10)/10
    swc.wilt = round(swc.q.con[[1]]*100)/100
    bucket.size=1000
  }else{
    df = gcc.met.pace.df
    swc.cap = 0.13
    swc.wilt = 0.05
    bucket.size=300
  }
  
  # par values####
  par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
    f.t.opt = c(10,25,40,NA,NA,NA),
    f.extract = c(0.2,1.5,8,NA,NA,NA),
    f.sec = c(0.01,0.15,0.5,NA,NA,NA),
    f.growth = c(0.01,0.15,0.3,NA,NA,NA),
    q = c(0.1,3,15,NA,NA,NA),
    q.s = c(0.1,1,15,NA,NA,NA))
  row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  
  # mcmc fitting
  fit.mcmc.2q.func(df,
                   n.iter = 50000,
                   species.in=species.vec[i],
                   prep.in = 'Control', temp.in ='Ambient',
                   my.fun = phenoGrass.func.v13,
                   out.nm.note='v1.2q.', 
                   use.smooth = TRUE,cal.initial = TRUE,day.lag = 3,
                   swc.capacity = swc.cap,swc.wilt = swc.wilt,
                   bucket.size = bucket.size,
                   par.df = par.df,q.given =NULL,q.s.given=NULL)
}

