# Fit model twith swc#################################
# this file reads data from PACE, DN, and Flux tower
# and fit CH model to them
#############################################################
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
  
  sd.swc <- sd(hufken.pace.pred$vwc,na.rm = T)
  
  resid.swc <- ((hufken.pace.pred$vwc.hufken - hufken.pace.pred$vwc)/sd.swc)^2
  
  return(c(sum(resid.gs,na.rm = T) + sum(resid.swc,na.rm = T)))
}
# ####
logLikelihood.func <- function (model.out){
  
  obs <- model.out$cover
  
  mod <- model.out$cover.hufken
  
  obs.swc <- model.out$vwc
  mod.swc <- model.out$vwc.hufken
  # if sd of obs can be estimated then use it 
  # otherwise assume a 10% sd
  if(sum(model.out$GCC.norm.sd,na.rm = T) != 0){
  obs.sd <- model.out$GCC.norm.sd 
  }else{
    obs.sd <- sd(model.out$GCC.norm,na.rm=T)
  }
  
  if(sum(model.out$vwc.sd,na.rm = T) != 0){
    obs.vwc.sd <- model.out$vwc.sd
  }else{
    obs.vwc.sd <- sd(model.out$vwc,na.rm=T)
  }
  
  
  logLi <-  - (0.5 * ((mod - obs)/obs.sd)^2 + log(obs.sd) + 0.5*log(2*pi))
  logLi.swc <- (0.5 * ((mod.swc - obs.swc)/obs.vwc.sd)^2 + log(obs.vwc.sd) + 0.5*log(2*pi))
  # mse <- mean((mod - obs)^2,na.rm=T)
  # logLi <- log(mse)
  return(sum(logLi,na.rm=T) + sum(logLi.swc,na.rm=T))
}
#read in data #################
ym.18.df <- get.ym.func(18)
ym.drt.df <- get.ym.func('Drought')
sd.df <- summaryBy(GCC ~ Date,
                   data = ym.drt.df,FUN=c(sd),na.rm=TRUE,keep.names = F)
ym.18.df$GCC.sd <- sd.df$GCC.sd
gcc.met.con.df <- get.paddock.func('control')

# loop through all spcies/site##################
species.vec <- c('Bis','Luc','Fes','Rye','ym','flux')
species.vec.nm <- c("Bis", "Med","Fes", "Lol", "YM" , "FT" )
for (i in seq_along(species.vec)){
  
  # use different soil water cap and wilt for different site
  if(species.vec[i]=='ym'){
    df = ym.18.df
    ym.met.df <- readRDS('cache/ym/ym.met.rds')
    
    swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    swc.cap = 0.3#round(swc.ym.con[[2]]*10)/10
    swc.wilt = .05#round(swc.ym.con[[1]]*100)/100
    bucket.size=1000
  }else if(species.vec[i]=='flux'){
    df = gcc.met.con.df
    swc.q.con <- quantile(gcc.met.con.df$vwc,na.rm=T,probs = c(0.01,0.99))
    swc.cap =  0.3#round(swc.q.con[[2]]*10)/10
    swc.wilt = 0.01#round(swc.q.con[[1]]*100)/100
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
                   n.iter = 1,
                   species.in=species.vec[i],
                   prep.in = 'Control', temp.in ='Ambient',
                   my.fun = phenoGrass.func.v13,
                   out.nm.note='.test.v1.2q.swc.', 
                   use.smooth = TRUE,cal.initial = TRUE,day.lag = 3,
                   swc.capacity = swc.cap,swc.wilt = swc.wilt,
                   bucket.size = bucket.size,
                   par.df = par.df,q.given =NULL,q.s.given=NULL)
}
#####
plot.ts.ci.swc.fit.func <- function(fn){
  
  # fn <-  paste0('tmp/pred.smv13.q1.qs0.chain.',species.vec[i],'.Control.Ambient.rds')
  # fn <- paste0('tmp/pred.smv13.q1.qs0.chain.Bis.Control.Ambient.rds')
  hufken.pace.pred <- readRDS(fn)
  
  hufken.pace.pred <- hufken.pace.pred[order(hufken.pace.pred$Date),]
  # hufken.pace.pred <- hufken.pace.pred[1:(nrow(hufken.pace.pred)),]
  
  # timeserie only
  
  
  if(sum(hufken.pace.pred$GCC.norm.sd,na.rm=T)!=0){
    hufken.pace.pred <- hufken.pace.pred[!is.na(hufken.pace.pred$GCC.norm.sd),]
    plot(GCC.norm~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
         xlab=' ',ylab='Cover',ylim=c(0,1),col = col.df$iris[4],
         xaxt='n')
    # add ci for obs
    hi.vec <- hufken.pace.pred$GCC.norm.smooth+hufken.pace.pred$GCC.norm.sd
    low.vec <- hufken.pace.pred$GCC.norm.smooth-hufken.pace.pred$GCC.norm.sd
    polygon(x = c(hufken.pace.pred$Date,
                  rev(hufken.pace.pred$Date)),
            y=c(hi.vec,rev(low.vec)),
            col=t_col(col.df$iris[4],80),border = NA
    )
  }else{
    plot(GCC.norm~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
         xlab=' ',ylab='Cover',ylim=c(0,1),col = col.df$iris[4],
         xaxt='n')
    
    sd.gcc <- 0.08#sd(hufken.pace.pred$GCC.norm,na.rm=T)
    hi.vec <- hufken.pace.pred$GCC.norm.smooth+sd.gcc
    low.vec <- hufken.pace.pred$GCC.norm.smooth-sd.gcc
    polygon(x = c(hufken.pace.pred$Date,
                  rev(hufken.pace.pred$Date)),
            y=c(hi.vec,rev(low.vec)),
            col=t_col(col.df$iris[4],80),border = NA)
  }
  
  
  
  # plot model pred
  points(cover.hufken~Date,data = hufken.pace.pred,type='l',lwd='2',col=col.df$auLandscape[2],lty='solid')
  
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


}
#####
plot.ts.swc.swc.fit.func <- function(fn){
  
  # fn <-  paste0('tmp/pred.smv13.q1.qs0.chain.',species.vec[i],'.Control.Ambient.rds')
  # fn <- paste0('tmp/pred.smv13.q1.qs0.chain.Bis.Control.Ambient.rds')
  hufken.pace.pred <- readRDS(fn)
  
  hufken.pace.pred <- hufken.pace.pred[order(hufken.pace.pred$Date),]
  # hufken.pace.pred <- hufken.pace.pred[1:(nrow(hufken.pace.pred)),]
  
  # timeserie only
  
  # 
  # if(sum(hufken.pace.pred$GCC.norm.sd,na.rm=T)!=0){
  #   hufken.pace.pred <- hufken.pace.pred[!is.na(hufken.pace.pred$GCC.norm.sd),]
  #   plot(vwc~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
  #        xlab=' ',ylab='VWC',ylim=c(0.01,0.35),col = col.df$iris[4],
  #        xaxt='n')
  #   # add ci for obs
  #   hi.vec <- hufken.pace.pred$GCC.norm.smooth+hufken.pace.pred$GCC.norm.sd
  #   low.vec <- hufken.pace.pred$GCC.norm.smooth-hufken.pace.pred$GCC.norm.sd
  #   polygon(x = c(hufken.pace.pred$Date,
  #                 rev(hufken.pace.pred$Date)),
  #           y=c(hi.vec,rev(low.vec)),
  #           col=t_col(col.df$iris[4],80),border = NA
  #   )
  # }else{
  if(length(grep('Luc',fn))>0){
    plot(vwc~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
         xlab=' ',ylab='SWC',col = col.df$iris[4],ylim=c(0.01,0.2),
         xaxt='n') 
  }else{
    
    plot(vwc~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
         xlab=' ',ylab='SWC',col = col.df$iris[4],#ylim=c(0.01,0.35),
         xaxt='n') 
  }
  
  
  # sd.gcc <- hufken.pace.pred$vwc.sd #0.08#sd(hufken.pace.pred$GCC.norm,na.rm=T)
  # hi.vec <- hufken.pace.pred$GCC.norm.smooth+sd.gcc
  # low.vec <- hufken.pace.pred$GCC.norm.smooth-sd.gcc
  # polygon(x = c(hufken.pace.pred$Date,
  #               rev(hufken.pace.pred$Date)),
  #         y=c(hi.vec,rev(low.vec)),
  #         col=t_col(col.df$iris[4],80),border = NA)
  # }
  
  
  
  # plot model pred
  points(vwc.hufken~Date,data = hufken.pace.pred,type='l',lwd='2',col=col.df$auLandscape[2],lty='solid')
  
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
  # legend('topleft',legend = sprintf('(%s) %s',letters[i],species.vec.nm[i]),
  #        bty='n')
  
}
##################################################################
# make plots######################################################
# use different soil water cap and wilt for different site########
##################################################################
for (i in seq_along(species.vec)){
if(species.vec[i]=='ym'){
  df = ym.18.df
  # 
  ym.met.df <- readRDS('cache/ym/ym.met.rds')
  # 
  swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
  swc.cap = round(swc.ym.con[[2]]*10)/10
  swc.wilt = round(swc.ym.con[[1]]*100)/100
  bucket.size=1000
}else if(species.vec[i]=='flux'){
  df = gcc.met.con.df
  swc.q.con <- quantile(gcc.met.con.df$vwc,na.rm=T,probs = c(0.01,0.99))
  swc.cap =  0.3#round(swc.q.con[[2]]*10)/10
  swc.wilt = 0.01#round(swc.q.con[[1]]*100)/100
  bucket.size=1000
}else{
  df = gcc.met.pace.df
  swc.cap = 0.13
  swc.wilt = 0.05
  bucket.size=300
}
# 
par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
  f.t.opt = c(10,25,40,NA,NA,NA),
  f.extract = c(0.2,1.5,8,NA,NA,NA),
  f.sec = c(0.01,0.15,0.5,NA,NA,NA),
  f.growth = c(0.01,0.15,0.5,NA,NA,NA),
  q = c(0.1,3,15,NA,NA,NA),
  q.s = c(0.1,1,15,NA,NA,NA))
row.names(par.df) <- c('min','initial','max','fit','stdv','prop')

# do ploting and prediction for v11
plot.mcmc.func.2q(df,species.vec[i],
                  prep.in='Control',temp.in='Ambient',
                  my.fun = phenoGrass.func.v13,
                  nm.note='.test.v1.2q.swc.',use.smooth = TRUE,day.lag = 3,
                  swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,
                  bucket.size = bucket.size)
# plot.title.func(species.vec[i])
}
# fn.v1.swc <- paste0('cache/sm.test.v1.2q.swc.chain.',species.vec[i],'.Control.Ambient.rds')
# df.v1.swc <- readRDS(fn.v1.swc)

###############
r.value.list <- list()
pdf('figures/withSWC.pdf',width = 5*2,height = 5*2*.618)
par(mar=c(5,5,1,1),mfrow=c(2,2))
# species.vec <- c('Bis','Luc','Fes','Rye','ym','flux')
species.vec <- c('ym','flux')
species.vec.nm <- c('YM','FT')
for (i in seq_along(species.vec)){
  if(species.vec[i]=='ym'){
    df = ym.18.df
    # 
    ym.met.df <- readRDS('cache/ym/ym.met.rds')
    # 
    swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    swc.cap = round(swc.ym.con[[2]]*10)/10
    swc.wilt = round(swc.ym.con[[1]]*100)/100
    bucket.size=1000
  }else if(species.vec[i]=='flux'){
    df = gcc.met.con.df
    swc.q.con <- quantile(gcc.met.con.df$vwc,na.rm=T,probs = c(0.01,0.99))
    swc.cap =  0.3#round(swc.q.con[[2]]*10)/10
    swc.wilt = 0.01#round(swc.q.con[[1]]*100)/100
    bucket.size=1000
  }else{
    df = gcc.met.pace.df
    swc.cap = 0.13
    swc.wilt = 0.05
    bucket.size=300
  }
  # 
  par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
    f.t.opt = c(10,25,40,NA,NA,NA),
    f.extract = c(0.2,1.5,8,NA,NA,NA),
    f.sec = c(0.01,0.15,0.5,NA,NA,NA),
    f.growth = c(0.01,0.15,0.5,NA,NA,NA),
    q = c(0.1,3,15,NA,NA,NA),
    q.s = c(0.1,1,15,NA,NA,NA))
  row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  fn <-  paste0('tmp/pred.sm.test.v1.2q.swc.chain.',species.vec[i],'.Control.Ambient.rds')
  df.swc.pred <- readRDS(fn)
  plot.ts.ci.swc.fit.func(fn)
  # 
  fn.v1 <- paste0('cache/smv1.2q.chain.',species.vec[i],'.Control.Ambient.rds')
  df.v1 <- readRDS(fn.v1)
  pred.df <- plot.with.given.par.func(df = df,fit.par.vec = df.v1[[1]][1,],
                                      species.in = species.vec[i], 
                                      prep.in='Control',temp.in='Ambient',
                                      my.fun = phenoGrass.func.v13,
                                      nm.note=' ',use.smooth = TRUE,day.lag = 3,
                                      swc.in.cap = swc.cap,
                                      swc.in.wilt = swc.wilt,
                                      bucket.size = bucket.size)
  points(cover.hufken~Date,data = pred.df,type='l',lty='dashed',lwd=3,col=col.df$auLandscape[2])
  if(i==1){
    legend('bottom',legend = c('With GCC','With GCC & SWC'),col= col.df$auLandscape[2],lty=c('dashed','solid'),bty='n',horiz = T)
  }
  # 
  legend('topleft',legend = sprintf('(%s) %s',letters[2*i-1],species.vec.nm[i]),
         bty='n')
  
  # r.value.list[[i]] <- data.frame(spc = species.vec.nm[i],
  #                          r.gcc =  get.r.func(pred.df$cover.hufken,
  #                                              pred.df$cover),
  #                          r.gcc.sec = get.r.func(df.swc.pred$cover.hufken,
  #                                       df.swc.pred$cover))
  # 
  plot.ts.swc.swc.fit.func(fn)
  points(vwc.hufken~Date,data = pred.df,type='l',lty='dashed',lwd=3,col=col.df$auLandscape[2])
  
  legend('topleft',legend = sprintf('(%s) %s',letters[2*i],species.vec.nm[i]),
         bty='n')
}
dev.off()

r.value.df <- do.call(rbind,r.value.list)

write.csv(r.value.df,'cache/r.swc.fit.csv',row.names = F)
