# plot.vec <- c(3,4,5,6)
# par(mfrow=c(3,1))
# for(plot.i in plot.vec){
#   dig.df <- dat.both.df
#   with(dig.df,plot(GCC.mean.con~Date,
#                    main = species.vec[plot.i],ylim=c(0,1)))
#   with(dig.df,points(GCC.mean.drt~Date,pch=16,col='red'))
#   
#   with(dig.df,points(cover.pred.con~Date,type='l',col='grey'))
#   with(dig.df,points(cover.pred.drt~Date,type='l',col='red'))
#   # par(new=T)
#   with(dig.df,plot(rain.con~Date,type='s'))
#   with(dig.df,plot(rain.drt~Date,type='s',col='red'))
#   
# }
pdf('figures/drt.ts.pdf',width = 8,height = 8*.618)
day.lag <- 3
source('r/pace_data_process.R')
source('r/ym_data_process.R')

# all data####
species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Rye','YM')
for (i in seq_along(species.vec)) {
  fn.con <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.Ambient.rds',
                    species.vec[i])
  # dat.con <- get.pace.func(gcc.met.pace.df,
  #                          species.in=species.vec[i],
  #                          prep.in = 'Control',
  #                          temp.in ='Ambient',
  #                          subplot = NA,
  #                          norm.min.max = NULL)
  
  dat.con <- readRDS(fn.con)
  dat.con <- dat.con[,c('Date','GCC.norm',
                        'GCC.norm.sd','swc.hufken','cover.hufken','Rain','harvest')]
  names(dat.con) <- c("Date","GCC.mean.con",
                      'GCC.sd.con',"vwc.mean.con",'cover.pred.con','rain.con','harvest')
  # get dry obs
  if( species.vec[i] != 'YM'){
    df <- get.pace.func(gcc.met.pace.df,
                        species.in=species.vec[i],
                        prep.in = 'Drought',
                        temp.in ='Ambient',
                        subplot = NA,
                        norm.min.max = NULL)
  }else{
    df <- get.pace.func(get.ym.func('Drought'),
                        species.in='ym',
                        prep.in = 'Drought',
                        temp.in ='Ambient',
                        subplot = NA,
                        norm.min.max = NULL)
  }
  dat.drought.obs <- df
  
  dat.drought.obs <- dat.drought.obs[,c('Date','GCC.norm',
                                        'GCC.norm.sd')]
  
  # 
  fn.drt <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.predict.Ambient.rds',
                    species.vec[i])
  
  dat.drought.pred <- readRDS(fn.drt)
  dat.drought.pred <- dat.drought.pred[,c('Date','swc.hufken','cover.hufken','Rain')]
  dat.drought <- merge(dat.drought.obs,dat.drought.pred,by='Date',all.y=T)
  
  names(dat.drought) <- c("Date","GCC.mean.drt",
                          'GCC.sd.drt',"vwc.mean.drt",'cover.pred.drt','rain.drt')
  
  dat.both.df <- merge(dat.con,dat.drought,
                       by = c('Date'))
# make plot
  dig.df <- dat.both.df
  with(dig.df,plot(GCC.mean.con~Date,
                   main = species.vec[plot.i],ylim=c(0,1),
                   xaxt='n',xlab=''))
  with(dig.df,points(GCC.mean.drt~Date,pch=16,col='red'))
  # 
  date.range = range(dig.df$Date,na.rm=T)
  mons.vec =  seq(date.range[1],date.range[2],by='mon')
  
  mon.c <- format(mons.vec,'%m')
  axis(1,at = mons.vec,labels = mon.c)
  # mtext('2018',side = 1,adj=0,line = 3)
  # mtext('2019',side = 1,adj=0.5,line = 3)
  yr.vec <- unique(year(dig.df$Date))
  where.c <-which(mon.c =='01') / length(mon.c)
  num.yr <- length(where.c)
  mtext(yr.vec[(length(yr.vec) - num.yr + 1):length(yr.vec)],side = 1,adj = where.c,line = 3)
  # 
  
  with(dig.df,points(cover.pred.con~Date,type='l',col='grey'))
  with(dig.df,points(cover.pred.drt~Date,type='l',col='red'))
  # # par(new=T)
  # with(dig.df,plot(rain.con~Date,type='s'))
  # with(dig.df,plot(rain.drt~Date,type='s',col='red'))
  
  # add harvest
  
  # clip(min(dig.df$Date), max(dig.df$Date), 0.0, 0.1)
  abline(v = dig.df$Date[dig.df$harvest ==1],lty='dotted')
  
}
dev.off()