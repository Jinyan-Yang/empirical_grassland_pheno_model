day.lag <- 3
source('r/pace_data_process.R')
source('r/ym_data_process.R')
source('r/read_spc_nm.R')
library(lubridate)

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
# 

# if(treat == 'Drought'){
#   df.14 <- readRDS('cache/ym/gcc_14.rds')
#   df.14$SubplotID <- 14
#   df.27 <- readRDS('cache/ym/gcc_27.rds')
#   df.27$SubplotID <- 27
#   df.38 <- readRDS('cache/ym/gcc_38.rds')
#   df.38$SubplotID <- 38
# plot(Rain_mm_Tot~Date,data = ym.14.df)
ym.14.df <- get.pace.func(get.ym.func(14),
                          species.in='ym',
                          prep.in = 'Drought',
                          temp.in ='Ambient',
                          subplot = NA,
                          norm.min.max = NULL)
ym.27.df <- get.pace.func(get.ym.func(27),
                          species.in='ym',
                          prep.in = 'Drought',
                          temp.in ='Ambient',
                          subplot = NA,
                          norm.min.max = NULL)
ym.38.df <- get.pace.func(get.ym.func(38),
                          species.in='ym',
                          prep.in = 'Drought',
                          temp.in ='Ambient',
                          subplot = NA,
                          norm.min.max = NULL)

ym.drt.df <- get.pace.func(get.ym.func('Drought'),
                          species.in='ym',
                          prep.in = 'Drought',
                          temp.in ='Ambient',
                          subplot = NA,
                          norm.min.max = NULL)
  # 

# ym.drought.df <- get.ym.func('Drought')
# ym.drought.df.sum <- get.pace.func(ym.drought.df,
#                                    species.in='ym',
#                                    prep.in = 'Drought',
#                                    temp.in ='Ambient',
#                                    subplot = NA,
#                                    norm.min.max = NULL)
# ym.drought.df.sum <- ym.drought.df.sum[,c('Date','GCC.norm',
#                                           'GCC.norm.sd','vwc')]
# names(ym.drought.df.sum) <- c("Date","GCC.mean.drt",
#                               'GCC.sd.drt',"vwc.mean.drt")
# 
# ym.both.df <- merge(ym.con.df.sum,ym.drought.df.sum,
#                     by = c('Date'))

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
big.rain.date <- as.Date('2020-2-8')
mean.con <- mean(ym.con.df.sum$GCC.mean.con[ym.con.df.sum$Date<big.rain.date],na.rm=T)
mean.annual.reduction <- c(mean(ym.14.df$GCC.norm[ym.14.df$Date<big.rain.date]) / mean.con,
                           mean(ym.27.df$GCC.norm[ym.27.df$Date<big.rain.date]) / mean.con,
                           mean(ym.38.df$GCC.norm[ym.38.df$Date<big.rain.date]) / mean.con)


mean.ym <- median(ym.drt.df$GCC.norm[ym.drt.df$Date<big.rain.date] / mean.con)
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
  
  dat.sum.wide$drt.iri <-  (dat.sum.wide$irrig.tot.Drought) /
    (dat.sum.wide$irrig.tot.Control)
  
  dat.sum.wide$species <- species.vec.nm[i] 
  dat.sum.wide$mean.response <- median(dat.sum.wide$drt.gcc,na.rm = T)
  pace.ls[[i]] <-  dat.sum.wide
  
}
pace.ls[[9]] <- data.frame(Shelter  = 1,
                           GCC.Control = NA,
                           irrig.tot.Control = NA,
                           GCC.Drought = NA,
                           irrig.tot.Drought = NA,
                           drt.gcc = mean.annual.reduction,
                           drt.iri = NA,
                           mean.response = mean.ym,
                           species = 'YM')
pace.effect.df <- do.call(rbind,pace.ls)
# 
pace.obs.mean.df <- pace.effect.df[,c('species','mean.response')]
pace.obs.mean.df <- pace.obs.mean.df[!duplicated(pace.obs.mean.df),]
pace.obs.mean.df <- pace.obs.mean.df[order(pace.obs.mean.df$mean.response),]
# 
pace.effect.df$species <- factor(pace.effect.df$species,
                                 levels = c(as.character(pace.obs.mean.df$species)))

# pace.effect.df$species <- factor(pace.effect.df$species,
#                                  levels = c(species.vec.nm[1:9]))


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
                                                mean(cover.pred.drt,na.rm=T) / 
                                                  mean(cover.pred.con,na.rm=T))
  # dat.both.df$mean.annual.reduction.rain <- with(dat.both.df,
  #                                                sum(rain.drt,na.rm=T) / 
  #                                                  sum(rain.con,na.rm=T))
  
  dat.both.df$spc <- species.vec.nm[i] 
  
  pace.model.ls[[i]] <-  dat.both.df
  
}

pace.model.df <- do.call(rbind,pace.model.ls)
pace.model.df$spc.factoir <- as.numeric(factor(pace.model.df$spc,levels = c(as.character(pace.obs.mean.df$species))))




# pace.model.df$spc.factoir <- droplevels(pace.model.df$spc.factoir)
# 
pdf('figures/plot.drt.shelter.pdf',width = 8,height = 8*.618)

# plot(c(drt.gcc)~species,data = pace.effect.df,xlab='',
#      ylab='Cover drought / Cover control',ylim=c(0,1),
#      col=c(1,1,2,2,2,3,3,3,4),
#      pch='')
plot(c(drt.gcc)~species,data = pace.effect.df,xlab='',
     ylab='Cover drought / Cover control',ylim=c(0.3,1),
     col=c(1,3,1,2,
           2,4,3,3,2),
     pch='')
# abline(h=1,lwd=2,col='grey',lty='dashed')
points(mean.annual.reduction.ped~spc.factoir,data = pace.model.df,
       pch=16,cex=2,col='grey30')
# plot(drt.iri~species,data = pace.effect.df)
# abline(h=1,lwd=2,col='grey',lty='dashed')

# #
# ci.df <- ci.df[ci.df$spc !='Flux Tower',]
# ci.df$plot.factor <- factor(ci.df$spc,levels = c(as.character(pace.obs.mean.df$species)))
# par(new=T)
# points(q.5 ~ plot.factor,data =ci.df)
# par(new=T)
# points(c(q.s/20)~plot.factor,data = par.val.df,ylim=c(0,20),pch=15,col='coral')
# points(c(q/20)~plot.factor,data = par.val.df,ylim=c(0,20),col='darkseagreen',pch=15)

par.val.df <- read.csv('cache/fittedParValue.csv')
par.val.df$spc <- species.vec.nm
par.val.df <- par.val.df[par.val.df$site != 'flux',]


par.val.df$plot.factor <- factor(par.val.df$spc,
                                            levels = pace.obs.mean.df$species)
# points(c(q.s/20)~plot.factor,data = par.val.df,ylim=c(0,20),pch=15,col='coral')
# points(c(q/20)~plot.factor,data = par.val.df,ylim=c(0,20),col='darkseagreen',pch=15)

dev.off()

