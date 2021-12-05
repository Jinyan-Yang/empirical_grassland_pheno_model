# Fit model to all data sets#################################
# this file reads data from PACE, DN, and Flux tower
# and fit CH model to them
#############################################################

# #this bit is not used but kept for debuging purpose
# source('r/load.R')
# source('r/read_spc_nm.R')

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
    f.t.opt = c(5,20,40,NA,NA,NA),
    f.extract = c(1,5,10,NA,NA,NA),
    f.sec = c(0.1,0.15,0.5,NA,NA,NA),
    f.growth = c(0.1,0.15,0.5,NA,NA,NA),
    q = c(0.1,3,15,NA,NA,NA),
    q.s = c(0.1,1,15,NA,NA,NA))
  row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  
  # # mcmc fitting
  # fit.mcmc.2q.func(df,
  #                  n.iter = 50000,
  #                  species.in=species.vec[i],
  #                  prep.in = 'Control', temp.in ='Ambient',
  #                  my.fun = phenoGrass.func.v13,
  #                  out.nm.note='smv13.2q.', 
  #                  use.smooth = TRUE,cal.initial = TRUE,day.lag = 3,
  #                  swc.capacity = swc.cap,swc.wilt = swc.wilt,
  #                  bucket.size = bucket.size,
  #                  par.df = par.df,q.given =NULL,q.s.given=NULL)
  # de opt fitting
  fit.mcmc.2q.func(df,
                   n.iter = 50000,
                   species.in=species.vec[i],
                   prep.in = 'Control', temp.in ='Ambient',
                   my.fun = phenoGrass.func.v13,
                   out.nm.note='de.v13.2q.',
                   use.smooth = TRUE,cal.initial = TRUE,day.lag = 3,
                   swc.capacity = swc.cap,swc.wilt = swc.wilt,
                   bucket.size = bucket.size,
                   par.df = par.df,q.given =NULL,q.s.given=NULL,use.mcmc=FALSE,de.note = '.new.t.')
  
  
  
}
pdf('figures/de.cover.pdf',width = 8,height = 8*.618)

for(spc.i in seq_along(species.vec)){
  
  if(species.vec[spc.i]=='ym'){
    df = ym.18.df
    ym.met.df <- readRDS('cache/ym/ym.met.rds')
    
    swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    swc.cap = round(swc.ym.con[[2]]*10)/10
    swc.wilt = round(swc.ym.con[[1]]*100)/100
    bucket.size=1000
  }else if(species.vec[spc.i]=='flux'){
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
  
  gcc.met.pace.df.16 <- get.pace.func(df,
                                      species.in = species.vec[spc.i],
                                      prep.in='Control',temp.in='Ambient',
                                      norm.min.max=NULL)
  # 
  par.fn <- sprintf('tmp/deopt_.new.t.%sControlAmbient.rds',species.vec[spc.i])
  de.par.vec <- readRDS(par.fn)
  # 
  hufken.pace.pred <- phenoGrass.func.v13(gcc.met.pace.df.16,
                                          f.h = 222,
                                          f.t.opt = de.par.vec[1],
                                          f.extract = de.par.vec[2],
                                          f.sec= de.par.vec[3],
                                          f.growth = de.par.vec[4],
                                          q =  de.par.vec[5],
                                          q.s =  de.par.vec[6],
                                          bucket.size = bucket.size,
                                          swc.wilt = swc.wilt ,
                                          swc.capacity = swc.cap ,
                                          t.max = 45,
                                          day.lay = 3,
                                          use.smooth = TRUE)
  
  # 
  # plot obs cover
  par(mar=c(5,5,1,5))
  plot(cover~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
       xlab=' ',ylab=expression(f[cover]),ylim=c(0,1),col = col.df$iris[4],
       xaxt='n')
  
  date.range = range(hufken.pace.pred$Date,na.rm=T)
  mons.vec =  seq(date.range[1],date.range[2],by='mon')
  
  mon.c <- format(mons.vec,'%m')
  axis(1,at = mons.vec,labels = mon.c)
  
  yr.vec <- unique(year(hufken.pace.pred$Date))
  where.c <-which(mon.c =='01') / length(mon.c)
  num.yr <- length(where.c)
  mtext(yr.vec[(length(yr.vec) - num.yr + 1):length(yr.vec)],side = 1,adj = where.c,line = 3)
  
  # plot model pred
  points(cover.hufken~Date,data = hufken.pace.pred,type='l',lwd='2',col=col.df$auLandscape[2],lty='solid')
  
  legend('topleft',legend = species.vec[spc.i],bty='n')
  
}

dev.off()