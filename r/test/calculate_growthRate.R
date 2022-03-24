# read data from Ak
fit.result.df <- readRDS('r/test/Polytunnel_experiment_jim/script/fit_gcc_swc.rds')
fit.result.df$Photo_Path <- as.factor(fit.result.df$Photo_Path)
# read fitted params from phenoMod
pars.df <- read.csv('cache/fittedParValue.csv')
pars.df <- pars.df[,c('site','f.growth','f.sec','q','q.s')]

# get fgrowth from harvest data
# DM = cover * 3629.82 +861.67)
# above relationship estimated from Amy's harvest data
# 
cInDryMass <- 0.5 #g/g
area.plot <- 4 #m2
pars.df$r.max.harvest <- (pars.df$f.growth * 3629.82 +861.67)*cInDryMass /area.plot 

# get fgrowth by assuming cover = fapar
pars.df$lai = log(1 - (pars.df$f.growth)) / (-0.5)
# get sla
sla.df <- read.csv('r/test/SLA_Nov2018.csv')
cm2m <- 1e-4
sla.df$lma <- sla.df$Target1leafDM / (sla.df$Target1LA_cm2*cm2m)#g/m2

# get spc mean
library(doBy)
sla.sum.df <- summaryBy(lma~SpeciesType,data = sla.df,FUN=mean,na.rm=T,
                    keep.names=T)
pars.sla.df <- merge(pars.df,sla.sum.df,
                     by.x = 'site',by.y = 'SpeciesType',all.x=T)

# 
pars.sla.df$r.max.fapar <- pars.sla.df$lai * pars.sla.df$lma *0.5
# 
pars.sla.df$pft <- c('leg','C4','C3','mix','C4','leg','C3','C4','C3','mix')
pars.sla.df$pft <- as.factor(pars.sla.df$pft)

mean(pars.sla.df$r.max.fapar,na.rm=T)
sd(pars.sla.df$r.max.fapar,na.rm=T)
# # 
# mean(pars.sla.df$r.max.fapar[pars.sla.df$pft=='c4'])
# sd(pars.sla.df$r.max.fapar[pars.sla.df$pft=='c4'])
# mean(pars.sla.df$r.max.fapar[pars.sla.df$pft=='c3'])
# sd(pars.sla.df$r.max.fapar[pars.sla.df$pft=='c3'])

# save data
write.csv(pars.sla.df,'r/test/growth.rate.csv',row.names = F,sep = ",")


# make plot to check data####
pdf('r/test/value compare.pdf',width = 8,height = 8)
par(mfrow=c(2,2),
    mar=c(3,5,1,1))
plot(f.sec~pft,data = pars.sla.df,xlab='',ylim=c(0.,0.5))
legend('topleft',legend='(a) PACE,YM & FT',bty='n')
# 
plot(r.s~Photo_Path,data = fit.result.df,
     xlab='',ylab='f.sec',ylim=c(0.,0.5))
legend('topleft',legend='(b) AK',bty='n')

# 
plot(q.s~pft,data = pars.sla.df,ylim=c(0,15),
     xlab='',ylab='q.sec')
legend('topleft',legend='(c) PACE,YM & FT',bty='n')
# 
plot(q.s~Photo_Path,data = fit.result.df,ylim=c(0,15),
     xlab='',ylab='q.sec')
legend('topleft',legend='(d) AK',bty='n')

# 
par(mfrow=c(2,1),
    mar=c(3,5,1,1))

plot(r.max.fapar~pft,data = pars.sla.df,
     xlab='',ylab='r.growth')
legend('topleft',legend='(a) ',bty='n')

# 
pars.sla.df$site <- as.factor(pars.sla.df$site)
plot(r.max.fapar~site,data = pars.sla.df,
     xlab='',ylab='r.growth')
legend('topleft',legend='(b) ',bty='n')

dev.off()


# 

# read LAi data from pace####
# this is not used for now
library(readxl)
lai.df <- read_xlsx('r/test/NonDestructiveSummary 2020-08-19 (1).xlsx', sheet = 'LAI')


lai.spc.df <- merge(lai.df,gcc.met.pace.df[,c('Date','Species','SubplotID',"Precipitation","Temperature",'GCC')],
                    by.x = c('Date','SubPlotID'),
                    by.y = c('Date','SubplotID'),all.x=T)

lai.spc.df$gcc.norm <- (lai.spc.df$GCC - 0.3) / 0.12
# 
lai.spc.df <- lai.spc.df[lai.spc.df$Precipitation == 'Control'&
                           lai.spc.df$Temperature == 'Ambient',]

plot(LAI~gcc.norm,data = lai.spc.df[lai.spc.df$Species == 'Kan',])
summary(lm(LAI~gcc.norm,data = lai.spc.df))

plot(LAI~Date,data = lai.spc.df[lai.spc.df$Species == 'Kan',],ylim=c(0,1))
points(gcc.norm~Date,data = lai.spc.df[lai.spc.df$Species == 'Kan',],pch=16,col='darkseagreen')
