######################################
# funtions to plot#
######################################
#
plot.title.func=function(species.in,where.to=0.5){
  par(mfrow=c(1,1),new=T,mar=c(1,1,1,1))
  plot(0,ann=F,axes=F,pch=' ')
  title(main = species.in,line = 0,adj=where.to)
}

# 
plot.mcmc.func.2q = function(df = gcc.met.pace.df,
                             species.in,prep.in,temp.in,subplot=NULL,
                             nm.note='',use.smooth=FALSE,
                             my.fun = phenoGrass.func.v11,
                             swc.in.cap = 0.13,swc.in.wilt = 0.05,bucket.size =300,
                             norm.min.max=NULL,
                             day.lag=3,
                             do.predict =NULL,
                             burn.proportion = 0.75,q.s.in=NULL,q.in=NULL){
  
  if(is.null(subplot)){
    gcc.met.pace.df.16 <- get.pace.func(df,
                                        species.in = species.in,
                                        prep.in = prep.in,
                                        temp.in =temp.in,
                                        norm.min.max=norm.min.max)

    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    if(is.null(do.predict)){
      fn=paste0('cache/',sm.nm,nm.note,'chain.',
                species.in,'.',prep.in,'.',temp.in,'.rds')
      rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',
                      species.in,'.',
                      prep.in,'.',
                      temp.in,'.rds')
      print(rds.nm)
    }else{
      fn=paste0('cache/',sm.nm,nm.note,'chain.',
                species.in,'.',do.predict,'.',temp.in,'.rds')
      rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',species.in,'.',do.predict,'.predict.',temp.in,'.rds')
      
      # reduce rainfall for drt species
      
      if(tolower(species.in) == 'ym'){
        gcc.met.pace.df.16$irrig.tot <- gcc.met.pace.df.16$irrig.tot * 0.35
        gcc.met.pace.df.16$Rain <- gcc.met.pace.df.16$Rain * 0.35
  
      }else{
        gcc.met.pace.df.16$irrig.tot[month(gcc.met.pace.df.16$Date) %in% 6:11] <-
          gcc.met.pace.df.16$irrig.tot[month(gcc.met.pace.df.16$Date) %in% 6:11] * 0.4
        gcc.met.pace.df.16$Rain[month(gcc.met.pace.df.16$Date) %in% 6:11] <-
          gcc.met.pace.df.16$Rain[month(gcc.met.pace.df.16$Date) %in% 6:11] * 0.4
      }
    }
    
   
    # fn=paste0('cache/chain.',species.in,'.',prep.in,'.',temp.in,'.rds')
    
    
  }else{
    species.in = subplot
    prep.in = ''
    temp.in =''
    gcc.met.pace.df.16 <- get.pace.func(gcc.met.pace.df,subplot = subplot)
    
    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    fn=paste0('cache/',sm.nm,nm.note,'chain.',subplot,'.rds')
    
    rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',subplot,'.rds')
  }
  report.df <<- gcc.met.pace.df.16
  # fn='cache/smv10.testchain.Fes.Control.Ambient.rds'
  # gcc.met.pace.df.16 <- gcc.met.pace.df.16[(gcc.met.pace.df.16$Date) < as.Date('2019-11-26'),]
  # gcc.met.pace.df.16$map <- 760
  
  # chain.fes <- readRDS('cache/chain.Rye.Control.Ambient.rds')
  # read chains 
  print(paste0('par file used: ',fn))
  in.chain =  readRDS(fn)
  # 
  # if(is.list(in.chain)){
  #   # assuming 1/3 burn in
  #   burnIn = 1
  #   chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round((burn.proportion)*nrow(m.in)):nrow(m.in),])
  #   
  #   chain.fes <- do.call(rbind,chain.3.ls.new)
  # }else{
  #   burnIn = nrow(in.chain)/3
  #   chain.fes <- in.chain
  # }

  # # # see how it works#####
  # par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
  #   f.t.opt = c(10,15,20,NA,NA,NA),
  #   f.extract = c(0.05,0.075,0.1,NA,NA,NA),
  #   f.sec = c(0.05,0.1,0.15,NA,NA,NA),
  #   f.growth = c(0.1,0.2,0.3,NA,NA,NA),
  #   q = c(0.001,1,2,NA,NA,NA),
  #   q.s = c(0.001,1,2,NA,NA,NA))
  # row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  
  # par.df["fit",] <- colMeans(chain.fes[burnIn:nrow(chain.fes),])
  
  # fit.par.vec <- apply(chain.fes[burnIn:nrow(chain.fes),],2,median)

  # # read tge best performing par
  # fn = paste0('tmp/bic.',sm.nm,nm.note,'chain.',
  #             species.in,'.',
  #             prep.in,'.',
  #             temp.in,'.rds')
  # in.chain =  readRDS(fn)
  # fit.par.vec <- in.chain[which(in.chain$bic==min(in.chain$bic)),]
  # 
  chain.fes <- do.call(rbind,in.chain)
  fit.par.vec <- subset(chain.fes[which(chain.fes$ll==max(chain.fes$ll)),],select=-c(ll))
  fit.par.vec <- fit.par.vec[!duplicated(fit.par.vec),]
  # 
  
  if(!'q' %in% names(fit.par.vec)){
    print(paste0('sensitivities of growth set to ',c(q.in)))
    fit.par.vec$q = q.in
  }
  if(!'q.s' %in% names(fit.par.vec)){
    print(paste0('sensitivities of senesence set to ',c(q.s.in)))
    fit.par.vec$q.s = q.s.in
  }
  
  # if(length(fit.par.vec)<5){
  #   fit.par.vec[5:6] <-c(q.in,q.s.in)
  #   print(paste0('sensitivities of growth and senesence set to ',c(q.in,q.s.in)))
  # }else if(length(fit.par.vec)<6){
  #   fit.par.vec[6] <- q.s.in
  #   print(paste0('sensitivitiy of senesence set to 1',q.s.in))
  # }
  

  # par.df["fit",] <- fit.par.vec
  print(fit.par.vec)
  # 
  # bucket.size = 300
  gcc.met.pace.df.16 <- gcc.met.pace.df.16[order(gcc.met.pace.df.16$Date),]
  hufken.pace.pred <- my.fun(gcc.met.pace.df.16,
                             f.h = 222,
                             f.t.opt = fit.par.vec$f.t.opt,
                             f.extract = fit.par.vec$f.extract,
                             f.sec= fit.par.vec$f.sec,
                             f.growth = fit.par.vec$f.growth,
                             q =  fit.par.vec$q,
                             q.s =  fit.par.vec$q.s,
                             bucket.size = bucket.size,
                             swc.wilt = swc.in.wilt ,
                             swc.capacity = swc.in.cap ,
                             t.max = 45,
                             day.lay = day.lag,use.smooth = use.smooth)
  
  # save prediction for future use
  saveRDS(hufken.pace.pred,rds.nm)

  # plot swc
  par(mfrow=c(2,2))
  par(mar=c(5,5,1,5))
  plot(vwc.hufken~Date,data = hufken.pace.pred,type='s',
       ann=F,axes=F,col = col.df$bushBySea[3],ylim=c(0,0.3),lwd='2')
  points(vwc~Date,data = hufken.pace.pred,type='s',lty='dashed',
         col = col.df$bushBySea[3],lwd='2')
  
  legend('topleft',legend = c('OBS','MOD'),lty = c('dashed','solid'),col='black')
  
  max.irrig = round(max(hufken.pace.pred$irrig.tot,na.rm=T))
  step.tmp <- floor((swc.in.cap /5)*100)/100
  axis(2,at = seq(0,0.3,by=step.tmp),labels = seq(0,0.3,by=step.tmp))
  mtext('VWC',side = 2,line = 3)
  
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

  # plot obs cover
  par(mar=c(5,5,1,5))
  plot(cover~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
       xlab=' ',ylab=expression(f[cover]),ylim=c(0,1),col = col.df$iris[4],
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
  points(cover.hufken~Date,data = hufken.pace.pred,type='l',lwd='2',col=col.df$auLandscape[2],lty='solid')
  # plot irrig
  par(new=TRUE)
  plot(irrig.tot~Date,data = hufken.pace.pred,type='s',
       ann=F,axes=F,col = 'navy')
  max.irrig = ceiling(max(hufken.pace.pred$irrig.tot,na.rm=T))
  axis(4,at = seq(0,max.irrig,by=10),labels = seq(0,max.irrig,by=10))
  mtext('irrigation (mm)',side = 4,line = 3)
  
  # legend('topleft',legend = paste0(species.in,prep.in,temp.in),bty='n')
  
  # clip(min(hufken.pace.pred$Date), max(hufken.pace.pred$Date), 0.0, 0.1)
  # abline(v = hufken.pace.pred$Date[hufken.pace.pred$harvest ==1],lty='dotted')
  
  # par(new=T)
  # 
  # plot(irrig.tot~Date,data = hufken.pace.pred,type='s',
  #      ann=F,axes=F,col = 'lightskyblue')
  # max.irrig = round(max(hufken.pace.pred$irrig.tot,na.rm=T))
  # axis(4,at = seq(0,max.irrig,by=10),labels = seq(0,max.irrig,by=10))
  # mtext('irrigation (mm)',side = 4)
  
  # vwc scater
  plot(vwc~vwc.hufken,data = hufken.pace.pred,pch=16,col='grey',
       xlab='MOD_VWC',ylab='OBS_VWC')
  abline(a=0,b=1)
  
  # scatter plot
  plot(cover~cover.hufken,data = hufken.pace.pred,pch=16,col='grey',
       xlab='MOD_GCC',ylab='OBS_GCC')
  abline(a=0,b=1)
  
}

# this will do for a modis site ######
# not used for now
plot.mcmc.func.2q.modis = function(df = gcc.met.pace.df,
                             species.in,prep.in,temp.in,subplot=NULL,
                             nm.note='',use.smooth=FALSE,
                             my.fun = phenoGrass.func.v11,
                             swc.in.cap = 0.13,swc.in.wilt = 0.05,bucket.size =300,
                             norm.min.max=NULL,time.series = TRUE,day.lag=1){
  
  if(is.null(subplot)){
    gcc.met.pace.df.16 <- get.pace.func(df,
                                        species.in = species.in,
                                        prep.in = prep.in,
                                        temp.in =temp.in,
                                        norm.min.max=norm.min.max)
    
    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    fn=paste0('cache/',sm.nm,nm.note,'chain.',
              species.in,'.',prep.in,'.',temp.in,'.rds')
    rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',species.in,'.',prep.in,'.',temp.in,'.rds')
    # fn=paste0('cache/chain.',species.in,'.',prep.in,'.',temp.in,'.rds')
    
    
  }else{
    species.in = subplot
    prep.in = ''
    temp.in =''
    gcc.met.pace.df.16 <- get.pace.func(gcc.met.pace.df,subplot = subplot)
    
    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    fn=paste0('cache/',sm.nm,nm.note,'chain.',subplot,'.rds')
    rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',subplot,'.rds')
  }
  
  # read chains 
  in.chain =  readRDS(fn)
  
  if(is.list(in.chain)){
    # assuming 1/3 burn in
    burnIn = 1
    chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round(3*nrow(m.in)/4):nrow(m.in),])
    
    chain.fes <- do.call(rbind,chain.3.ls.new)
  }else{
    burnIn = nrow(in.chain)/3
    chain.fes <- in.chain
  }

  # # see how it works#####
  # par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
  #   f.t.opt = c(10,15,20,NA,NA,NA),
  #   f.extract = c(0.05,0.075,0.1,NA,NA,NA),
  #   f.sec = c(0.05,0.1,0.15,NA,NA,NA),
  #   f.growth = c(0.1,0.2,0.3,NA,NA,NA),
  #   q = c(0.001,1,2,NA,NA,NA),
  #   q.s = c(0.001,1,2,NA,NA,NA))
  # row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  # par.df["fit",] <- colMeans(chain.fes[burnIn:nrow(chain.fes),])
  # par.df["fit",] <- colMeans(luc.d.a.df[burnIn:nrow(luc.d.a.df),])
  # 
  # bucket.size = 300
  hufken.pace.pred <- my.fun(gcc.met.pace.df.16,
                             f.h = 222,
                             f.t.opt = par.df["fit",1],
                             f.extract = par.df["fit",2],
                             f.sec= par.df["fit",3],
                             f.growth = par.df["fit",4],
                             q =  par.df["fit",5],
                             q.s =  par.df["fit",6],
                             bucket.size = bucket.size,
                             swc.wilt = swc.in.wilt ,
                             swc.capacity = swc.in.cap ,
                             t.max = 45,
                             day.lay = day.lag,use.smooth = use.smooth)
  
  # save prediction for future use
  saveRDS(hufken.pace.pred,rds.nm)
  
  # if(time.series==TRUE){
    # par(mar=c(5,5,1,5))
    par(mfrow=c(2,1))
    # plot obs cover
    par(mar=c(5,5,1,5))
    plot(cover~Date,data = hufken.pace.pred,type='p',pch=16,#lwd='2',
         xlab=' ',ylab=expression(f[cover]),ylim=c(0,1),col = col.df$iris[4],
         xaxt='n')
    
    # add dates
    date.range = range(hufken.pace.pred$Date,na.rm=T)
    mons.vec =  seq(date.range[1],date.range[2],by='mon')
    
    mon.c <- format(mons.vec,'%m')
    axis(1,at = mons.vec,labels = mon.c)

    yr.vec <- unique(year(hufken.pace.pred$Date))
    where.c <-which(mon.c =='01') / length(mon.c)
    num.yr <- length(where.c)
    mtext(yr.vec[(length(yr.vec) - num.yr + 1):length(yr.vec)],side = 1,adj = where.c,line = 3)
    
    # plot model pred
    points(cover.hufken~Date,data = hufken.pace.pred,type='l',lwd='2',col=col.df$auLandscape[2],lty='dashed')
    
    # plot rainfall
    par(new=T)
    plot(irrig.tot~Date,data = hufken.pace.pred,type='s',
         ann=F,axes=F,col = 'lightskyblue')
    max.irrig = round(max(hufken.pace.pred$irrig.tot,na.rm=T))
    axis(4,at = seq(0,max.irrig,by=10),labels = seq(0,max.irrig,by=10))
    mtext('irrigation (mm)',side = 4,line=3)
  # }else{
    # scatter plot
    plot(cover~cover.hufken,data = hufken.pace.pred,pch=16,col='grey',
         xlab='MOD_GCC',ylab='OBS_GCC')
    abline(a=0,b=1)
  # }
  
}

# functions to ckeck if fitting makes sense##############
plot.check.mcmc.func=function(chain.in,species.in='',nm.vec = c('Topt','f.extract','senescence','growth','q','qs'),
                              plot.par.no=1:6,col.in='grey'){
  
  burnIn = round(nrow(chain.in) / 4*3)
  
  par(mfrow=c(length(plot.par.no)/2,2),mar=c(5,5,1,1))
  
  for(i in plot.par.no){
    hist(chain.in[burnIn:nrow(chain.in),i],xlab = nm.vec[i],
         main='',col=col.in)
    
  }
  # plot.title.func(species.in = species.in)
}

# 
# functions to chist##############
plot.hist.func=function(chain.in,
                        nm.vec = c('Topt','f.extract','senescence','growth','q','qs'),
                        plot.par.no=1:6){
  
  burnIn = round(nrow(chain.in[[1]]) / 4*3)
  
  par(mfrow=c(length(plot.par.no)/2,2),mar=c(5,5,1,1))

  range.list <- list(c(10,40),
                     c(0.2,8),
                     c(0.01,0.5),
                     c(0.01,0.5),
                     c(0.01,15),
                     c(0.01,15))
  for(i in plot.par.no){
    x.range <- range.list[[i]]
    hist(chain.in[[1]][burnIn:nrow(chain.in[[1]]),i],xlab = nm.vec[i],
         main='',col=t_col(palette()[1],percent = 70),border =NA,xlim = x.range)
    hist(chain.in[[2]][burnIn:nrow(chain.in[[2]]),i],xlab = nm.vec[i],
         main='',col=t_col(palette()[2],percent = 70),add=TRUE,border =NA)
    hist(chain.in[[3]][burnIn:nrow(chain.in[[3]]),i],xlab = nm.vec[i],
         main='',col=t_col(palette()[3],percent = 70),add=TRUE,border =NA)
    
  }
  # plot.title.func(species.in = species.in)
}

# 
acceptance.func <- function(vec){
  sum(!duplicated(vec)) / length(vec)
}

# 
plot.line.mcmc.func <- function(chain.3.ls,val.nm,
                                nm.vec = c('Topt','f.extract','senescence','growth','q','qs'),
                                range.iter = NULL){
  
  # 
  n.iter <- nrow(chain.3.ls[[1]])
  
  if(is.null(range.iter)){
    range.iter <- round(0.5*n.iter):n.iter
  }
  
  # 
  min.val <- min(min(chain.3.ls[[1]][range.iter,val.nm]),min(chain.3.ls[[3]][range.iter,val.nm]),min(chain.3.ls[[2]][range.iter,val.nm]))
  
  max.val <- max(max(chain.3.ls[[1]][range.iter,val.nm]),max(chain.3.ls[[3]][range.iter,val.nm]),max(chain.3.ls[[2]][range.iter,val.nm]))
  # 
  plot(chain.3.ls[[1]][range.iter,val.nm],pch=16,ylim=c(min.val,max.val),ylab=nm.vec[val.nm],xlab='Iteration',col=1)
  print(acceptance.func(chain.3.ls[[1]][range.iter,val.nm]))
  for (i in 2:3) {
    points(chain.3.ls[[i]][range.iter,val.nm],pch=16,col=i)
  }
}

# make colot transprent##############
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

# function to plot ts####

plot.ts.ci.func <- function(fn){
  
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

  
  if(sum(hufken.pace.pred$harvest,na.rm=T)>0){
    # add harvest
    clip(min(hufken.pace.pred$Date), max(hufken.pace.pred$Date), 0.0, 0.1)
    abline(v = hufken.pace.pred$Date[hufken.pace.pred$harvest ==1],lty='dotted')
  }



}

# ####
plot.with.given.par.func = function(df = gcc.met.pace.df,
                                    fit.par.vec,
                                    species.in,prep.in,temp.in,subplot=NULL,
                                    nm.note='',use.smooth=FALSE,
                                    my.fun = phenoGrass.func.v11,
                                    swc.in.cap = 0.13,swc.in.wilt = 0.05,bucket.size =300,
                                    norm.min.max=NULL,
                                    day.lag=3,q.s.in=NULL,q.in=NULL){
  do.predict=NULL
  if(is.null(subplot)){
    gcc.met.pace.df.16 <- get.pace.func(df,
                                        species.in = species.in,
                                        prep.in = prep.in,
                                        temp.in =temp.in,
                                        norm.min.max=norm.min.max)
    
    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    if(is.null(do.predict)){
      fn=paste0('cache/',sm.nm,nm.note,'chain.',
                species.in,'.',prep.in,'.',temp.in,'.rds')
      rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',
                      species.in,'.',
                      prep.in,'.',
                      temp.in,'.rds')
      print(rds.nm)
    }else{
      fn=paste0('cache/',sm.nm,nm.note,'chain.',
                species.in,'.',do.predict,'.',temp.in,'.rds')
      rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',species.in,'.',do.predict,'.predict.',temp.in,'.rds')
      
      # reduce rainfall for drt species
      
      if(tolower(species.in) == 'ym'){
        gcc.met.pace.df.16$irrig.tot <- gcc.met.pace.df.16$irrig.tot * 0.35
        gcc.met.pace.df.16$Rain <- gcc.met.pace.df.16$Rain * 0.35
        
      }else{
        gcc.met.pace.df.16$irrig.tot[month(gcc.met.pace.df.16$Date) %in% 6:11] <-
          gcc.met.pace.df.16$irrig.tot[month(gcc.met.pace.df.16$Date) %in% 6:11] * 0.4
        gcc.met.pace.df.16$Rain[month(gcc.met.pace.df.16$Date) %in% 6:11] <-
          gcc.met.pace.df.16$Rain[month(gcc.met.pace.df.16$Date) %in% 6:11] * 0.4
      }
    }
    
    
    # fn=paste0('cache/chain.',species.in,'.',prep.in,'.',temp.in,'.rds')
    
    
  }else{
    species.in = subplot
    prep.in = ''
    temp.in =''
    gcc.met.pace.df.16 <- get.pace.func(gcc.met.pace.df,subplot = subplot)
    
    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    fn=paste0('cache/',sm.nm,nm.note,'chain.',subplot,'.rds')
    
    rds.nm = paste0('tmp/pred.',sm.nm,nm.note,'chain.',subplot,'.rds')
  }
  report.df <<- gcc.met.pace.df.16
  # fn='cache/smv10.testchain.Fes.Control.Ambient.rds'
  # gcc.met.pace.df.16 <- gcc.met.pace.df.16[(gcc.met.pace.df.16$Date) < as.Date('2019-11-26'),]
  # gcc.met.pace.df.16$map <- 760
  
  # chain.fes <- readRDS('cache/chain.Rye.Control.Ambient.rds')
  # # read chains 
  # print(paste0('par file used: ',fn))
  # in.chain =  readRDS(fn)
  # 
  # if(is.list(in.chain)){
  #   # assuming 1/3 burn in
  #   burnIn = 1
  #   chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round((burn.proportion)*nrow(m.in)):nrow(m.in),])
  #   
  #   chain.fes <- do.call(rbind,chain.3.ls.new)
  # }else{
  #   burnIn = nrow(in.chain)/3
  #   chain.fes <- in.chain
  # }
  
  # # # see how it works#####
  # par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
  #   f.t.opt = c(10,15,20,NA,NA,NA),
  #   f.extract = c(0.05,0.075,0.1,NA,NA,NA),
  #   f.sec = c(0.05,0.1,0.15,NA,NA,NA),
  #   f.growth = c(0.1,0.2,0.3,NA,NA,NA),
  #   q = c(0.001,1,2,NA,NA,NA),
  #   q.s = c(0.001,1,2,NA,NA,NA))
  # row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  
  # par.df["fit",] <- colMeans(chain.fes[burnIn:nrow(chain.fes),])
  
  # fit.par.vec <- apply(chain.fes[burnIn:nrow(chain.fes),],2,median)
  
  # # read tge best performing par
  # fn = paste0('tmp/bic.',sm.nm,nm.note,'chain.',
  #             species.in,'.',
  #             prep.in,'.',
  #             temp.in,'.rds')
  # in.chain =  readRDS(fn)
  # fit.par.vec <- in.chain[which(in.chain$bic==min(in.chain$bic)),]
  # # 
  # chain.fes <- do.call(rbind,in.chain)
  # fit.par.vec <- subset(chain.fes[which(chain.fes$ll==max(chain.fes$ll)),],select=-c(ll))
  # fit.par.vec <- fit.par.vec[!duplicated(fit.par.vec),]
  # 
  
  if(!'q' %in% names(fit.par.vec)){
    print(paste0('sensitivities of growth set to ',c(q.in)))
    fit.par.vec$q = q.in
  }
  if(!'q.s' %in% names(fit.par.vec)){
    print(paste0('sensitivities of senesence set to ',c(q.s.in)))
    fit.par.vec$q.s = q.s.in
  }
  
  # if(length(fit.par.vec)<5){
  #   fit.par.vec[5:6] <-c(q.in,q.s.in)
  #   print(paste0('sensitivities of growth and senesence set to ',c(q.in,q.s.in)))
  # }else if(length(fit.par.vec)<6){
  #   fit.par.vec[6] <- q.s.in
  #   print(paste0('sensitivitiy of senesence set to 1',q.s.in))
  # }
  
  
  # par.df["fit",] <- fit.par.vec
  print(fit.par.vec)
  # 
  # bucket.size = 300
  gcc.met.pace.df.16 <- gcc.met.pace.df.16[order(gcc.met.pace.df.16$Date),]
  hufken.pace.pred <- my.fun(gcc.met.pace.df.16,
                             f.h = 222,
                             f.t.opt = fit.par.vec$f.t.opt,
                             f.extract = fit.par.vec$f.extract,
                             f.sec= fit.par.vec$f.sec,
                             f.growth = fit.par.vec$f.growth,
                             q =  fit.par.vec$q,
                             q.s =  fit.par.vec$q.s,
                             bucket.size = bucket.size,
                             swc.wilt = swc.in.wilt ,
                             swc.capacity = swc.in.cap ,
                             t.max = 45,
                             day.lay = day.lag,use.smooth = use.smooth)
  
  return(hufken.pace.pred)
  
}
