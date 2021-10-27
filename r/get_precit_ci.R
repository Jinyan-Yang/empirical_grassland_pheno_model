
day.lag <- 3
source('r/function_get_ci.R')
source('r/pace_data_process.R')
source('r/ym_data_process.R')
source('r/v13_common_fun.R')
source('models/hufkens/hufkensV13.R')
source('r/process_paddock_gcc_met.R')
# devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
library(zoo)
# source('r/plot.mcmc.r')

ym.18.df <- get.ym.func(18)
gcc.met.con.df <- get.paddock.func('control')
# species.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','ym','flux')
species.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','ym','flux')

for (i in seq_along(species.vec)) {
  
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
  
  get.mod.ci.func(df,species.vec[i],'Control','Ambient',
                  my.fun = phenoGrass.func.v13,
                  nm.note='smv13.2q.',use.smooth = TRUE,day.lag = 3,
                  swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,
                  bucket.size = bucket.size,
                  sample.size = 200)

}
# # do the same for v10
# for (i in seq_along(species.vec)) {
#   
#   # use different soil water cap and wilt for different site
#   if(species.vec[i]=='ym'){
#     df = ym.18.df
#     # 
#     c.wd <- getwd()
#     setwd('c:/repo/dn_gcc/')
#     ym.met.df <- readRDS('cache/ym.met.rds')
#     setwd(c.wd)
#     # 
#     swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
#     swc.cap = round(swc.ym.con[[2]]*10)/10
#     swc.wilt = round(swc.ym.con[[1]]*100)/100
#     bucket.size=1000
#   }else if(species.vec[i]=='flux'){
#     df = gcc.met.con.df
#     swc.q.con <- quantile(gcc.met.con.df$vwc,na.rm=T,probs = c(0.01,0.99))
#     swc.cap =  round(swc.q.con[[2]]*10)/10
#     swc.wilt = round(swc.q.con[[1]]*100)/100
#     bucket.size=1000
#   }else{
#     df = gcc.met.pace.df
#     swc.cap = 0.13
#     swc.wilt = 0.05
#     bucket.size=300
#   }
#   
#   get.mod.ci.func(df,species.vec[i],'Control','Ambient',
#                   my.fun = phenoGrass.func.v13,
#                   nm.note='v13.q1.qs0.',use.smooth = TRUE,day.lag = 3,
#                   swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,
#                   bucket.size = bucket.size,q.s.in=0,q.in=1)
#   
# }


