plot.ts.func <- function(hufken.pace.pred){
  # plot irrig
  
  plot(irrig.tot~Date,data = hufken.pace.pred,type='s',
       ann=F,axes=F,col = 'navy')
  max.irrig = ceiling(max(hufken.pace.pred$irrig.tot,na.rm=T))
  axis(4,at = seq(0,max.irrig,by=10),labels = seq(0,max.irrig,by=10))
  mtext('irrigation (mm)',side = 4,line = 3)
  par(new=TRUE)
  plot(GCC.norm~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
       xlab=' ',ylab='Cover',ylim=c(0,1),col = col.df$iris[4],
       xaxt='n')
  
  
  date.range = range(hufken.pace.pred$Date,na.rm=T)
  mons.vec =  seq(date.range[1],date.range[2],by='mon')
  
  mon.c <- format(mons.vec,'%m')
  axis(1,at = mons.vec,labels = mon.c)
  # mtext('2018',side = 1,adj=0,line = 3)
  # mtext('2019',side = 1,adj=0.5,line = 3)
  yr.vec <- unique(year(hufken.pace.pred$Date))
  where.c <-which(mon.c =='01') / length(mon.c)
  num.yr <- length(where.c)
  mtext(yr.vec[(length(yr.vec) - num.yr + 1):length(yr.vec)],side = 1,adj = where.c,line = 3)
  
  # plot model pred
  points(cover.50~Date,data = hufken.pace.pred,type='l',lwd='2',col=col.df$auLandscape[2],lty='solid')
  
}

ym.18.df <- get.ym.func(18)

gcc.met.con.df <- get.paddock.func('control')

# pdf('figures/obs_pred.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec)) {
  
  # use different soil water cap and wilt for different site
  if(species.vec[i]=='ym'){
    df = ym.18.df
    # 
    # c.wd <- getwd()
    # setwd('c:/repo/dn_gcc/')
    ym.met.df <- readRDS('cache/ym/ym.met.rds')
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
  # do ploting and prediction for v11
  
  gcc.met.pace.df.16 <- get.pace.func(df,
                                      species.in = species.vec[i],
                                      prep.in = 'Control',
                                      temp.in ='Ambient',
                                      norm.min.max=NULL)
  
  
  # 
  gcc.met.pace.df.16 <- gcc.met.pace.df.16[order(gcc.met.pace.df.16$Date),]
  
  # read de param
  fn <- sprintf('tmp/deopt_%sControlAmbient.rds', species.vec[i])
  param.vec <- readRDS(fn)
  # 
  hufken.de.pred <- phenoGrass.func.v13(gcc.met.pace.df.16,
                             f.h = 222,
                             f.t.opt = param.vec[1],
                             f.extract = param.vec[2],
                             f.sec= param.vec[3],
                             f.growth = param.vec[4],
                             q =  param.vec[5],
                             q.s =  param.vec[6],
                             bucket.size = bucket.size,
                             swc.wilt = swc.wilt ,
                             swc.capacity = swc.cap ,
                             t.max = 45,
                             day.lay = 3,use.smooth = TRUE)
  
# 
  # par(mar=c(5,5,1,5))
  fn <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  hufken.pace.pred <- readRDS(fn)
  
  ci.fm <- sprintf('tmp/ci.smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  ci.m <- readRDS(ci.fm)
  hufken.pace.pred$cover.05 <- ci.m[1,]
  hufken.pace.pred$cover.95 <- ci.m[2,]
  hufken.pace.pred$cover.50 <- ci.m[3,]
  
  plot.ts.func(hufken.pace.pred)
  points(cover.hufken~Date,hufken.de.pred,type='l',lty='dotted',col='red')
  # 
  legend('topleft',legend = species.vec.nm[i])

}
