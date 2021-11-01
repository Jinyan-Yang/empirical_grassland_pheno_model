# set up the environment####
# this lag day term is outdated but need to be here to use the older version of code
day.lag <- 3
source('r/pace_data_process.R')
source('r/ym_data_process.R')
source('r/v13_common_fun.R')
source('models/hufkens/hufkensV13.R')
source('r/process_paddock_gcc_met.R')

library(zoo)
library(foreach)
library(doParallel)

#read in data 
ym.18.df <- get.ym.func(18)
gcc.met.con.df <- get.paddock.func('control')

species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','ym','flux')
# species.vec <- c('Kan','ym','flux')
# species.vec <- c('Bis','Luc','Dig','Rho','Fes','Pha','Rye','flux')
# loop through all spcies/site
for (i in seq_along(species.vec)){
  
  # use different soil water cap and wilt for different site
  if(species.vec[i]=='ym'){
    df = ym.18.df
    # 
    c.wd <- getwd()
    setwd('c:/repo/dn_gcc/')
    ym.met.df <- readRDS('cache/ym.met.rds')
    setwd(c.wd)
    # 
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
    f.t.opt = c(5,20,40,NA,NA,NA),
    f.extract = c(0.01,0.5,0.8,NA,NA,NA),
    f.sec = c(0.001,0.05,0.5,NA,NA,NA),
    f.growth = c(0.001,0.01,0.5,NA,NA,NA))
  row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  
  # mcmc fitting
  fit.mcmc.2q.func(df,
                   n.iter = 50000,
                   species.in=species.vec[i],
                   prep.in = 'Control', temp.in ='Ambient',
                   my.fun = phenoGrass.func.v10,
                   out.nm.note='v10.', 
                   use.smooth = TRUE,cal.initial = TRUE,day.lag = 3,
                   swc.capacity = swc.cap,swc.wilt = swc.wilt,bucket.size = bucket.size,
                   par.df = par.df,q.given =1,q.s.given=0)
  
  
}
