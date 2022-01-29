source('r/plot.mcmc.r')
source('r/read_spc_nm.R')
# 
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
library(doBy)
library(lubridate)

# ts with ci####
pdf('figures/v10_ts.pdf',width = 5*2,height = 5*5*.618)
par(mfrow=c(5,2))
par(mar=c(5,5,1,5))
for (i in seq_along(species.vec)) {
  fn <-  paste0('tmp/pred.smv0.chain.',species.vec[i],'.Control.Ambient.rds')
  plot.ts.ci.func(fn)
  # # fn <- paste0('tmp/pred.smv13.q1.qs0.chain.Bis.Control.Ambient.rds')
  # hufken.pace.pred <- readRDS(fn)
  # # timeserie only
  # 
  # plot(GCC.norm~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
  #      xlab=' ',ylab='Cover',ylim=c(0,1),col = col.df$iris[4],
  #      xaxt='n')
  # 
  # # add date
  # date.range = range(hufken.pace.pred$Date,na.rm=T)
  # s.date <- as.Date(paste0(format(date.range[1],'%Y-%m'),'-01'))
  # e.date <- as.Date( paste0(format(date.range[2],'%Y-%m'),'-01'))
  # 
  # mons.vec =  seq(s.date,e.date,by='mon')
  # 
  # mon.c <- format(mons.vec,'%m')
  # axis(1,at = mons.vec,labels = mon.c)
  # # mtext('2018',side = 1,adj=0,line = 3)
  # # mtext('2019',side = 1,adj=0.5,line = 3)
  # yr.vec <- unique(year(hufken.pace.pred$Date))
  # where.c <- which(mon.c =='01') / length(mon.c)
  # num.yr <- length(where.c)
  # mtext(yr.vec[(length(yr.vec) - num.yr + 1):length(yr.vec)],side = 1,adj = where.c,line = 3)
  # 
  # 
  # # add ci
  # # ci.smv13.2q.07072021.chain.flux.Control.Ambient.rds
  # ci.fm <- sprintf('tmp/ci.smv13.q1.qs0.chain.%s.Control.Ambient.rds',
  #                  species.vec[i])
  # ci.m <- readRDS(ci.fm)
  # hufken.pace.pred$cover.05 <- ci.m[1,]
  # hufken.pace.pred$cover.95 <- ci.m[2,]
  # 
  # # plot model pred
  # points(ci.m[3,]~hufken.pace.pred$Date,type='l',lwd='2',col=col.df$auLandscape[2],lty='solid')
  # 
  # polygon(x = c(hufken.pace.pred$Date,
  #               rev(hufken.pace.pred$Date)),
  #         y=c(hufken.pace.pred$cover.95,rev(hufken.pace.pred$cover.05)),
  #         col=t_col(col.df$iris[4],60),border = NA
  # )
  # 
  # legend('topleft',legend = sprintf('(%s) %s',letters[i],species.vec.nm[i]),
  #        bty='n')
  
}

dev.off()
