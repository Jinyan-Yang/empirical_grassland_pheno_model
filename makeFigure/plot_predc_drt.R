day.lag <- 3
source('r/pace_data_process.R')
source('r/ym_data_process.R')
source('r/v13_common_fun.R')
source('models/hufkens/hufkensV13.R')
# source('r/process_paddock_gcc_met.R')
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
library(zoo)
library(lubridate)
source('r/plot.mcmc.r')
source('r/read_spc_nm.R')

ym.18.df <- get.ym.func('Drought')
# gcc.met.con.df <- get.paddock.func('control')
# species.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','ym')
# species.vec <- c('Bis','Luc',
#                  'Dig','Kan','Rho','Fes',
#                  'Rye','ym')
# species.vec <- c('Bis','Luc','Dig','Wal','Rho','Fes','Pha','Rye','ym')
# species.vec <- 'Kan'

# 
pdf('figures/pred_drought.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec[1:9])) {
  
  # use different soil water cap and wilt for different site
  if(species.vec[i]=='ym'){
    df = ym.18.df
    # df$irrig.tot <-  df$irrig.tot * 0.35
    # 
    # c.wd <- getwd()
    # setwd('c:/repo/dn_gcc/')
    ym.met.df <- readRDS('cache/ym/ym.met.rds')
    # setwd(c.wd)
    # 
    swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    swc.cap = round(swc.ym.con[[2]]*10)/10
    swc.wilt = round(swc.ym.con[[1]]*100)/100
    bucket.size=1000
  }else if(species.vec[i]=='flux'){
    stop('flux site should not be here')
  }else{
    df = gcc.met.pace.df
    # df$irrig.tot[month(df$Date) %in% 6:11] <- df$irrig.tot[month(df$Date) %in% 6:11] * 0.4
    # df$Rain[month(df$Date) %in% 6:11] <- df$Rain[month(df$Date) %in% 6:11] * 0.4
    # df$irrig.tot <- df$irrig.tot * 0.4
    swc.cap = 0.13
    swc.wilt = 0.05
    bucket.size=300
  }
  
  # plot control under control
  plot.mcmc.func.2q(df,species.vec[i],prep.in='Drought',temp.in='Ambient',
                    my.fun = phenoGrass.func.v13,
                    nm.note='smv13.2q.',use.smooth = TRUE,day.lag = 3,
                    swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,bucket.size = bucket.size,
                    do.predict = 'Control')
  
  plot.title.func(species.vec[i])
}
dev.off()

