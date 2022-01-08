############################################
# plot drought response and comare to model predictions
############################################

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

# filter data to before flood
big.rain.date <- as.Date('2020-2-8')
mean.con <- mean(ym.con.df.sum$GCC.mean.con[ym.con.df.sum$Date<big.rain.date],na.rm=T)
mean.annual.reduction <- c(mean(ym.14.df$GCC.norm[ym.14.df$Date<big.rain.date]) / mean.con,
                           mean(ym.27.df$GCC.norm[ym.27.df$Date<big.rain.date]) / mean.con,
                           mean(ym.38.df$GCC.norm[ym.38.df$Date<big.rain.date]) / mean.con)


mean.ym <- median(ym.drt.df$GCC.norm[ym.drt.df$Date<big.rain.date] / mean.con)

# pace data
pace.ls <- list()
# use only the late half of treatment
pace.sub.df <- gcc.met.pace.df[month(gcc.met.pace.df$Date) %in% 8:11,]
pace.sub.df <- pace.sub.df[pace.sub.df$Temperature == 'Ambient',]
for (i in seq_along(species.vec[1:8])) {

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


# get modelled responses
pace.model.ls <- list()
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
  
 
  dat.both.df$mean.annual.reduction.ped <- with(dat.both.df,
                                                mean(cover.pred.drt,na.rm=T) / 
                                                  mean(cover.pred.con,na.rm=T))
  
  dat.both.df$spc <- species.vec.nm[i] 
  
  pace.model.ls[[i]] <-  dat.both.df
  
}

pace.model.df <- do.call(rbind,pace.model.ls)
pace.model.df$spc.factoir <- as.numeric(factor(pace.model.df$spc,
                                               levels = c(as.character(pace.obs.mean.df$species))))

# plot 
pdf('figures/plot.drt.shelter.pdf',width = 8,height = 8*.618)
palette(c(col.df$iris))
# obs
plot(c(drt.gcc)~species,data = pace.effect.df,xlab='',
     ylab='Cover drought / Cover control',ylim=c(0.3,1),
     col=c(1,3,1,2,
           2,4,3,3,2),
     pch='')
# mod
points(mean.annual.reduction.ped~spc.factoir,data = pace.model.df,
       pch=16,cex=2,col='grey30')


dev.off()

