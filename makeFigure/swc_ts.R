source('r/plot.mcmc.r')
source('r/read_spc_nm.R')
# 
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
library(doBy)
library(lubridate)
swc.df <- readRDS('cache/pace/swc.rds')
treat.info.df <- gcc.met.pace.df[,c('SubplotID' ,'Temperature', 'Precipitation')]
treat.info.df <- treat.info.df[!duplicated(treat.info.df),]
swc.trt.df <- merge(swc.df,treat.info.df,by.x = 'SubPlotID',
                    by.y='SubplotID',all.x=T)
swc.trt.df <- swc.trt.df[swc.trt.df$Temperature =='Ambient' &
                           swc.trt.df$Precipitation == 'Control'&
                           swc.trt.df$Location=='Upper',]

# 
plot.swc.ts <- function(whereToLook){
  for (i in seq_along(species.vec)) {
    # '
    fn <-  paste0(whereToLook,species.vec[i],'.Control.Ambient.rds')
    # plot.ts.ci.func(fn)
    hufken.pace.pred <- readRDS(fn)
    # 
    hufken.pace.pred <- hufken.pace.pred[order(hufken.pace.pred$Date),]
    
    # 
    if(species.vec[i]=='ym'){
      swc.trt.df.spc <- readRDS('cache/ym/ym.met.rds')
      palette((col.df$beauty[c(5,3,1)]))
      plot(swc.5~Date,data = swc.trt.df.spc,pch=16,col = 1,
           xlab=' ',ylab='Soil volumetric water content',ylim=c(0,0.3),
           xaxt='n',cex=0.4)
      # points(swc.5~Date,data = swc.trt.df.spc,pch=16,col = 2,cex=0.2)
      points(swc.15.45~Date,data = swc.trt.df.spc,pch=16,col =2,cex=0.4)
      points(swc.45.75~Date,data = swc.trt.df.spc,pch=16,col = 3,cex=0.4)

      points(vwc.hufken~Date,data = hufken.pace.pred,type='l',lwd=4,col = col.df$auLandscape[2])
      legend('bottomright',legend = c('MOD','<5cm','15-45cm','45-75cm'),lty=c('solid',NA,NA,NA),pch=c(NA,16,16,16),
             bty='n',col = c(col.df$auLandscape[2],palette()),
             lwd=2,horiz = T)
      # plot(vwc.hufken~Date,data = hufken.pace.pred,
      #      type='l',lwd=4,col = col.df$auLandscape[2],
      #      xlab=' ',ylab='Soil volumetric water content',ylim=c(0,0.4),
      #      xaxt='n',cex=1)
      # 
      # points(vwc~Date,data = hufken.pace.pred,pch=16,col = t_col(col.df$iris[4],50))
      # legend('topright',legend = c('MOD','OBS'),lty=c('solid','dashed'),
      #        bty='n',col = c(col.df$auLandscape[2],t_col(col.df$iris[4],50)),
      #        lwd=2)
    }else if(species.vec[i]=='flux'){
      
      plot(vwc.hufken~Date,data = hufken.pace.pred,
           type='l',lwd=4,col = col.df$auLandscape[2],
           xlab=' ',ylab='Soil volumetric water content',ylim=c(0,0.4),
           xaxt='n',cex=1)
      
      points(vwc~Date,data = hufken.pace.pred,pch=16,col = t_col(col.df$iris[4],50))
      legend('topright',legend = c('MOD','OBS'),lty=c('solid','dashed'),
             bty='n',col = c(col.df$auLandscape[2],t_col(col.df$iris[4],50)),
             lwd=2)
    }else{
      swc.trt.df.spc <- swc.trt.df[swc.trt.df$Species == species.vec[i],]
      
      plot(vwc.hufken~Date,data = hufken.pace.pred,
           type='l',lwd=4,col = col.df$auLandscape[2],
           xlab=' ',ylab='Soil volumetric water content',ylim=c(0,0.2),
           xaxt='n',cex=1)
      if(nrow(swc.trt.df.spc)>0){
        swc.sum <- summaryBy(vwc ~ Date +Species,data = swc.trt.df.spc,
                             FUN=c(mean,sd),na.rm=T)
        hi.vec <- swc.sum$vwc.mean+swc.sum$vwc.sd
        low.vec <- swc.sum$vwc.mean-swc.sum$vwc.sd
        points(vwc.mean~Date,data = swc.sum,type='l',lty='dashed',lwd=3,col = t_col(col.df$iris[4],50))
        polygon(x = c(swc.sum$Date,
                      rev(swc.sum$Date)),
                y=c(hi.vec,rev(low.vec)),
                col=t_col(col.df$iris[4],80),border = NA)
        legend('topright',legend = c('MOD','OBS'),lty=c('solid','dashed'),
               bty='n',col = c(col.df$auLandscape[2],t_col(col.df$iris[4],50)),
               lwd=2)
        
      }
      legend('topright',legend = c('MOD'),lty=c('solid'),
             bty='n',col = col.df$auLandscape[2],
             lwd=2)
    }
    # add date
    date.range = range(hufken.pace.pred$Date,na.rm=T)
    s.date <- as.Date(paste0(format(date.range[1],'%Y-%m'),'-01'))
    e.date <- as.Date( paste0(format(date.range[2],'%Y-%m'),'-01'))
    
    mons.vec =  seq(s.date,e.date,by='mon')
    
    mon.c <- format(mons.vec,'%m')
    axis(1,at = mons.vec,labels = mon.c)
    # mtext('2018',side = 1,adj=0,line = 3)
    # mtext('2019',side = 1,adj=0.5,line = 3)
    yr.vec <- unique(year(hufken.pace.pred$Date))
    where.c <- which(mon.c =='01') / length(mon.c)
    num.yr <- length(where.c)
    mtext(yr.vec[(length(yr.vec) - num.yr + 1):length(yr.vec)],side = 1,adj = where.c,line = 3)
    
    # 
    legend('topleft',legend = sprintf('(%s) %s',letters[i],species.vec.nm[i]),
           bty='n')
    
    
    
  }
}

levels(swc.trt.df$Species) <- c("Bis","Fes","Luc","Rye")
# ts with ci####
pdf('figures/v10_swc_ts.pdf',width = 5*2,height = 5*5*.618)
par(mfrow=c(5,2))
par(mar=c(5,5,1,5))

plot.swc.ts(whereToLook = 'tmp/pred.smv0.chain.')

dev.off()

# ts with ci####
pdf('figures/v11_swc_ts.pdf',width = 5*2,height = 5*5*.618)
par(mfrow=c(5,2))
par(mar=c(5,5,1,5))

plot.swc.ts(whereToLook = 'tmp/pred.smv1.2q.chain.')
dev.off()