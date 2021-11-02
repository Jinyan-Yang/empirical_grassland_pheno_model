source('r/plot.mcmc.r')
source('r/read_spc_nm.R')
# 
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
library(doBy)
library(lubridate)

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

# species.vec <- c("Bis",    "Dig",  "Fes",    "Kan",    
#                  "Luc",  "Rho",    "Rye",
#                  'ym')
# species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM','flux')
# species.vec.plot.nm <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM','Flux Tower')
# #########################################
palette(c(col.df$iris,col.df$daisy))
png('figures/obs_fit_TS_scatter.png',height = 400*2,width = 400/.618)
par(mfrow =c(2,1))

# plot obs cover
par(mar=c(5,5,1,5))
fn <- 'tmp/pred.smsmv13.2q.chain.ym.Control.Ambient.rds'
hufken.pace.pred <- readRDS(fn)

ci.fm <- ('tmp/ci.smsmv13.2q.chain.ym.Control.Ambient.rds')
ci.m <- readRDS(ci.fm)
hufken.pace.pred$cover.05 <- ci.m[1,]
hufken.pace.pred$cover.95 <- ci.m[2,]
hufken.pace.pred$cover.50 <- ci.m[3,]

plot.ts.func(hufken.pace.pred)

polygon(x = c(hufken.pace.pred$Date,
              rev(hufken.pace.pred$Date)),
        y=c(hufken.pace.pred$cover.95,rev(hufken.pace.pred$cover.05)),
        col=t_col(col.df$iris[4],60),border = NA
        )
# 
# points(cover.05~Date,data = hufken.pace.pred,
#        type='l',lwd=2,col=col.df$iris[5],lty='dashed')
# points(cover.95~Date,data = hufken.pace.pred,
#        type='l',lwd=2,col= col.df$iris[5],lty='dashed')

legend('topleft',legend = '(a) YM',bty='n')
legend('topright',legend = c('OBS','MOD'),
       pch=c(16,NA),lty=c(NA,'solid'),
       col=c( col.df$iris[4],col.df$auLandscape[2]),
       bty='n')

for (i in seq_along(species.vec)){
  fn <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  hufken.pace.pred <- readRDS(fn)
  # 
  ci.fm <- sprintf('tmp/ci.smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  ci.m <- readRDS(ci.fm)
  # hufken.pace.pred$cover.05 <- ci.m[1,]
  # hufken.pace.pred$cover.95 <- ci.m[2,]
  hufken.pace.pred$cover.50 <- ci.m[3,]
  # 
  if(i == 1){
    plot(GCC.norm~cover.50,data = hufken.pace.pred,
         xlim=c(0,1),ylim=c(0,1),
         xlab='Modelled cover',ylab = 'Observed cover',pch=16,col=i)
   
  }else{
    points(GCC.norm~cover.50,data = hufken.pace.pred,
         xlim=c(0,1),ylim=c(0,1),
         pch=16,col=i)
  }
  legend('bottomright',legend = species.vec.nm,col=palette(),
         pch=16,bty='n')
  legend('topleft',legend = '(b)',bty='n')
  abline(a=0,b=1,lty='dashed',col='grey',lwd=2)

}

dev.off()

# ###########
# species.vec <- c("Bis","Dig", "Fes", "Kan",    
#                  "Luc",  "Rho",    "Rye")

# # species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye')
# species.vec <- plot.nm.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','Pha','YM')
# png('figures/fit_pred_TS_scatter.png',height = 400*2,width = 400/.618)
# par(mfrow =c(2,1))
# palette(c(col.df$iris,col.df$daisy))
# # plot obs cover
# par(mar=c(5,5,1,5))
# fn <- 'tmp/pred.smsmv13.2q.chain.ym.Control.predict.Ambient.rds'
# hufken.pace.pred <- readRDS(fn)
# 
# plot.ts.func(hufken.pace.pred)
# legend('topleft',legend = '(a) YM',bty='n')
# legend('topright',legend = c('OBS','MOD'),
#        pch=c(16,NA),lty=c(NA,'solid'),
#        col=c( col.df$iris[4],col.df$auLandscape[2]),
#        bty='n')
# 
# 
# for (i in seq_along(species.vec)){
#   fn <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.predict.Ambient.rds',species.vec[i])
#   hufken.pace.pred <- readRDS(fn)
# 
#   if(i == 1){
#     plot(GCC.norm~cover.hufken,data = hufken.pace.pred,
#          xlim=c(0,1),ylim=c(0,1),
#          xlab='MOD_cover',ylab = 'OBS_cover',pch=16,col=i)
#   }else{
#     points(GCC.norm~cover.hufken,data = hufken.pace.pred,
#            xlim=c(0,1),ylim=c(0,1),pch=16,col=i)
#   }
#   legend('bottomright',legend = plot.nm.vec,col=palette(),
#          pch=16,bty='n')
#   legend('topleft',legend = '(b)',bty='n')
#   abline(a=0,b=1,lty='dashed',col='grey',lwd=2)
# 
# }
# 
# dev.off()
# 
