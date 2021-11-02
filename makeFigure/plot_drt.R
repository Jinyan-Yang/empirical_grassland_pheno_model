day.lag <- 3
source('r/pace_data_process.R')
source('r/ym_data_process.R')
source('r/read_spc_nm.R')

devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
palette(c(col.df$iris))


# process ym data
ym.18.df <- get.ym.func(18)
ym.con.df.sum <- get.pace.func(ym.18.df,
                               species.in='ym',
                               prep.in = 'Control',
                               temp.in ='Ambient',
                               subplot = NA,
                               norm.min.max = NULL)
ym.con.df.sum <- ym.con.df.sum[,c('Date','GCC.norm',
                                  'GCC.norm.sd','vwc')]

names(ym.con.df.sum) <- c("Date","GCC.mean.con",
                          'GCC.sd.con',"vwc.mean.con")


ym.drought.df <- get.ym.func('Drought')
ym.drought.df.sum <- get.pace.func(ym.drought.df,
                                   species.in='ym',
                                   prep.in = 'Drought',
                                   temp.in ='Ambient',
                                   subplot = NA,
                                   norm.min.max = NULL)
ym.drought.df.sum <- ym.drought.df.sum[,c('Date','GCC.norm',
                                          'GCC.norm.sd','vwc')]
names(ym.drought.df.sum) <- c("Date","GCC.mean.drt",
                              'GCC.sd.drt',"vwc.mean.drt")

ym.both.df <- merge(ym.con.df.sum,ym.drought.df.sum,
                    by = c('Date'))
library(lubridate)
# ym.both.df <- ym.both.df[month(ym.both.df$Date) %in% 6:11,]
# 
# plot(GCC.mean.drt~Date,data = ym.both.df,type='l',col='red')
# points(GCC.mean.con~Date,data = ym.both.df,type='l',col='navy')
# # 
# plot(vwc.mean.drt~Date,data = ym.both.df,type='l',col='red')
# points(vwc.mean.con~Date,data = ym.both.df,type='l',col='navy')
# # 
# daily.reduction <- with(ym.both.df,GCC.mean.drt/GCC.mean.con)
# hist(daily.reduction,freq = F)
# 
mean.annual.reduction <- with(ym.both.df,(mean(GCC.mean.drt) ) / 
                                (mean(GCC.mean.con)))


# mean.annual.reduction.vwc <- with(ym.both.df,sum(vwc.mean.drt) / sum(vwc.mean.con))
# pace data
pace.ls <- list()
# species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye')

pace.sub.df <- gcc.met.pace.df[month(gcc.met.pace.df$Date) %in% 8:11,]
pace.sub.df <- pace.sub.df[pace.sub.df$Temperature == 'Ambient',]
for (i in seq_along(species.vec[1:8])) {
  # fn.con <- sprintf('tmp/pred.smv13.2qchain.%s.Control.Ambient.rds',
  #                   species.vec[i])
  # dat.con <- get.pace.func(gcc.met.pace.df,
  #                          species.in=species.vec[i],
  #                          prep.in = 'Control',
  #                          temp.in ='Ambient',
  #                          subplot = NA,
  #                          norm.min.max = NULL)
  
  dat.con <- pace.sub.df[pace.sub.df$Species == species.vec[i],]
  
  dat.sum <- summaryBy(GCC + irrig.tot~Precipitation + Shelter,data = dat.con,
                       FUN=mean,na.rm=T,keep.names =T)
  dat.sum <- dat.sum[complete.cases(dat.sum),]
  
  dat.sum.wide <- reshape(dat.sum, idvar = "Shelter", timevar = "Precipitation", direction = "wide")
  
  dat.sum.wide$drt.gcc <-  (dat.sum.wide$GCC.Drought -0.3) /
    (dat.sum.wide$GCC.Control -0.3)
  
  dat.sum.wide$drt.iri <-  (dat.sum.wide$irrig.tot.Drought-0.3) /
    (dat.sum.wide$irrig.tot.Control -0.3)
  
  dat.sum.wide$species <- species.vec.nm[i] 
  pace.ls[[i]] <-  dat.sum.wide
  
}
pace.ls[[9]] <- data.frame(Shelter  = 1,
                           GCC.Control = NA,
                           irrig.tot.Control = NA,
                           GCC.Drought = NA,
                           irrig.tot.Drought = NA,
                           drt.gcc = mean.annual.reduction,
                           drt.iri = NA,
                           species = 'YM')
pace.effect.df <- do.call(rbind,pace.ls)
pace.effect.df$species <- factor(pace.effect.df$species,
                                 levels = c(species.vec.nm))

# 
pace.model.ls <- list()
# species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM')
for (i in seq_along(species.vec[1:9])) {
  fn.con <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.Ambient.rds',
                    species.vec[i])
  
  dat.con <- readRDS(fn.con)
  dat.con <- dat.con[,c('Date','GCC.norm',
                        'GCC.norm.sd','vwc','cover.hufken','Rain')]
  names(dat.con) <- c("Date","GCC.mean.con",
                      'GCC.sd.con',"vwc.mean.con",'cover.pred.con','rain.con')
  
  fn.drt <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.predict.Ambient.rds',
                    species.vec[i])
  dat.drought <- readRDS(fn.drt)
  dat.drought <- dat.drought[,c('Date','GCC.norm',
                                'GCC.norm.sd','vwc','cover.hufken','Rain')]
  names(dat.drought) <- c("Date","GCC.mean.drt",
                          'GCC.sd.drt',"vwc.mean.drt",'cover.pred.drt','rain.drt')
  
  dat.both.df <- merge(dat.con,dat.drought,
                       by = c('Date'))
  if( species.vec[i] != 'YM'){
    dat.both.df <- dat.both.df[month(dat.both.df$Date) %in% 8:11,]
  }
  
  # # dat.both.df$daily.reduction.obs <- with(dat.both.df,
  #                                         GCC.mean.drt/GCC.mean.con)
  # dat.both.df$daily.reduction.pred <- with(dat.both.df,
  #                                          cover.pred.drt/cover.pred.con)
  # 
  # dat.both.df$mean.annual.reduction.obs <- with(dat.both.df,
  #                                               sum(GCC.mean.drt,na.rm=T) /
  #                                                 sum(GCC.mean.con,na.rm=T))
  dat.both.df$mean.annual.reduction.ped <- with(dat.both.df,
                                                sum(cover.pred.drt,na.rm=T) / 
                                                  sum(cover.pred.con,na.rm=T))
  # dat.both.df$mean.annual.reduction.rain <- with(dat.both.df,
  #                                                sum(rain.drt,na.rm=T) / 
  #                                                  sum(rain.con,na.rm=T))
  
  dat.both.df$spc <- species.vec.nm[i] 
  
  pace.model.ls[[i]] <-  dat.both.df
  
}

pace.model.df <- do.call(rbind,pace.model.ls)
pace.model.df$spc.factoir <- as.numeric(factor(pace.model.df$spc,levels = c(species.vec.nm)))
# 
pdf('figures/plot.drt.shelter.pdf',width = 8,height = 8*.618)

plot(c(drt.gcc)~species,data = pace.effect.df,xlab='',
     ylab='Reduction of cover under drought',ylim=c(0,1),
     col=c(1,1,2,2,2,3,3,3,4),
     pch='')
# abline(h=1,lwd=2,col='grey',lty='dashed')
points(mean.annual.reduction.ped~spc.factoir,data = pace.model.df,
       pch=16,cex=2,col='grey30')
# plot(drt.iri~species,data = pace.effect.df)
# abline(h=1,lwd=2,col='grey',lty='dashed')


dev.off()
