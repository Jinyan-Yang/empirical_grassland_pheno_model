# function to get the intersect of growth vs sensescence
solve.intersect.func <- function(vwc.in){
  swc.norm <-  (vwc.in- swc.wilt)/ (swc.capacity - swc.wilt)
  loss.f <- swc.norm^q
  loss.f.s <- (1-swc.norm)^q.s
  
  growth.vec <- 
    loss.f *
    (1 - cover.in / cover.max)
  
  senescence.vec <- (loss.f.s) *
    # (1 - cover.pred.vec[nm.day-1])*
    cover.in / cover.max
  
  return(growth.vec - senescence.vec)
}
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")

pdf('figures/theory_demo.pdf',width = 8,height = 2*8*.62)
par(mar=c(5,5,1,1))
par(mfrow=c(2,1))
# a########
# x <- seq(0,20,by=0.1)
# 
# plot.cal.function <- function(x,a,b=0.4){
#   y <- a*(exp(-b*x)-exp(-a*x)) / (a -b)
# }
# 
# 
# with(plot.df,
#      plot(y~x,data = plot.df[plot.df$x<=x[which.max(y)],],xlab='Time',ylab='Cover',type='l',col='black',lwd=3,
#           xaxt='n',yaxt='n',xlim=c(0,20)))
# 
# # plot(y~x,data = plot.df[plot.df$x<=x[which.max(y)],],xlab='Time',ylab='Cover',type='l',col='black',lwd=3,
# #      xaxt='n',yaxt='n',xlim=c(0,20))
# with(plot.df,
#      points(y~x,data = plot.df[plot.df$x>x[which.max(y)],],type='l',lty='dashed',lwd=3)
# )
# 
# with(plot.df,
#      points(y.3~x,data = plot.df[plot.df$x<=x[which.max(y.3)],],type='l',lwd=3,col='grey')
# )
# 
# with(plot.df,
#      points(y.3~x,data = plot.df[plot.df$x>x[which.max(y.3)],],type='l',lty='dashed',lwd=3,col='grey')
# )
# 
# # with(plot.df,
# #      points(y.4~x,data = plot.df[plot.df$x>x[which.max(y.4)],],type='l',lty='dashed',lwd=3,col='grey')
# # )
# # plot(y~x,xlab='Time',ylab='Cover',type='l',col='black',lwd=3,
# #      xaxt='n',yaxt='n',xlim=c(-5,25))
# # points(y.3~x,type='l',lty='dashed',lwd=3)
# 
# with(plot.df,
#      points(x[which.max(y)],y[which.max(y)],cex=3,pch=3))
# 
# with(plot.df,
#      points(x[which.max(y.3)],y.3[which.max(y.3)],cex=3,pch=3,col='grey'))
# # legend(-5,0.3,legend = 'High green-up sensitivity',bty='n')
# # legend(-5,0.05,legend = 'Low green-up sensitivity',bty='n')
# # 
# # legend(5,0.3,legend = 'High brown-down sensitivity',bty='n')
# # legend(5,0.05,legend = 'Low green-up sensitivity',bty='n')
# legend('topright',legend = c('Onset of Brown-down',
#                              'High green-up sensitivity', 
#                              'Low green-up sensitivity',
#                              'High brown-down sensitivity',
#                              'Low brown-down sensitivity'),bty='n',
#        pch=c(3,rep(NA,4)),
#        lty = c(NA,'solid','solid','dashed','dashed'),
#        col = c('black','black','grey','black','grey'))
source('models/hufkens/hufkensV13.R')
source('r/v13_common_fun.R')
day.lag <- 3
source('r/pace_data_process.R')
# source('r/ym_data_process.R')
# ym.18.df <- get.ym.func(18)
ym.processed.df <- get.pace.func(gcc.met.pace.df,
                                 species.in = 'Fes',
                                 prep.in = 'Control',
                                 temp.in ='Ambient',
                                 norm.min.max=NULL)

ym.processed.df <- ym.processed.df[year(ym.processed.df$Date)==2018,]

ym.processed.df$Tair <- 25
ym.processed.df$Tmax <- 30
ym.processed.df$Tmin <- 20
ym.processed.df$vwc <- 0.05
ym.processed.df$GCC.norm.smooth <- 0.0
# ym.processed.df$Rain_mm_Tot <- 0
# ym.processed.df$irrig.tot <- 0
# ym.processed.df$irrig.tot[10] <- 30
ym.processed.df$Rain <- 0
ym.processed.df$Rain[5] <- 50
# # 
# fn <- sprintf('cache/smv13.2q.chain.%s.Control.Ambient.rds','ym')
# in.chain =  readRDS(fn)
# burn.proportion=0.75
# if(is.list(in.chain)){
#   # assuming 1/3 burn in
#   burnIn = 1
#   chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round((1-burn.proportion)*nrow(m.in)):nrow(m.in),])
#   
#   chain.fes <- do.call(rbind,chain.3.ls.new)
# }else{
#   burnIn = nrow(in.chain)/3
#   chain.fes <- in.chain
# }
# # see how it works#####
par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
  f.t.opt = c(10,15,20,NA,NA,NA),
  f.extract = c(0.05,0.075,0.1,NA,NA,NA),
  f.sec = c(0.05,0.1,0.15,NA,NA,NA),
  f.growth = c(0.1,0.2,0.3,NA,NA,NA),
  q = c(0.001,1,2,NA,NA,NA),
  q.s = c(0.001,1,2,NA,NA,NA))
row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
# par.df["fit",] <- colMeans(chain.fes[burnIn:nrow(chain.fes),])

# par.df["fit",] <-apply(chain.fes[burnIn:nrow(chain.fes),],2,median)
par.df["fit",] <- c(26.01437,0.1029175,
                    0.1,0.1,
                    0.3143503,2.290821)
# 

hufken.pace.pred.lgls <- phenoGrass.func.v13(ym.processed.df,
                                             f.h = 222,
                                             f.t.opt = par.df["fit",1],
                                             f.extract = par.df["fit",2],
                                             f.sec= par.df["fit",3],
                                             f.growth = par.df["fit",4],
                                             q =  0.6,
                                             q.s =  4,
                                             bucket.size = 1000,
                                             swc.wilt = 0.05 ,
                                             swc.capacity = 0.3 ,
                                             t.max = 45,
                                             day.lay = 3,use.smooth = TRUE)


hufken.pace.pred.lghs <- phenoGrass.func.v13(ym.processed.df,
                                             f.h = 222,
                                             f.t.opt = par.df["fit",1],
                                             f.extract = par.df["fit",2],
                                             f.sec= par.df["fit",3],
                                             f.growth = par.df["fit",4],
                                             q =  0.6,
                                             q.s =  0.1,
                                             bucket.size = 1000,
                                             swc.wilt = 0.05 ,
                                             swc.capacity = 0.3 ,
                                             t.max = 45,
                                             day.lay = 3,use.smooth = TRUE)

hufken.pace.pred.hgls <- phenoGrass.func.v13(ym.processed.df,
                                             f.h = 222,
                                             f.t.opt = par.df["fit",1],
                                             f.extract = par.df["fit",2],
                                             f.sec= par.df["fit",3],
                                             f.growth = par.df["fit",4],
                                             q =  1,
                                             q.s =  4,
                                             bucket.size = 1000,
                                             swc.wilt = 0.05 ,
                                             swc.capacity = 0.3 ,
                                             t.max = 45,
                                             day.lay = 3,use.smooth = TRUE)

hufken.pace.pred.hghs <- phenoGrass.func.v13(ym.processed.df,
                                             f.h = 222,
                                             f.t.opt = par.df["fit",1],
                                             f.extract = par.df["fit",2],
                                             f.sec= par.df["fit",3],
                                             f.growth = par.df["fit",4],
                                             q =  1,
                                             q.s =  0.1,
                                             bucket.size = 1000,
                                             swc.wilt = 0.05 ,
                                             swc.capacity = 0.3 ,
                                             t.max = 45,
                                             day.lay = 3,use.smooth = TRUE)



palette(col.df$reptile)
palette(col.df$auLandscape[c(1,2,4,5)])
with(hufken.pace.pred.lgls,
     plot(cover.hufken~Date,type='l',lwd=3,col=1,
          xaxt='n',yaxt='n',xlab='Time',ylab = 'Cover'))


with(hufken.pace.pred.lghs,
     points(cover.hufken~Date,type='l',lwd=3,col=2))

with(hufken.pace.pred.hgls,
     points(cover.hufken~Date,type='l',lwd=3,col=3))

with(hufken.pace.pred.hghs,
     points(cover.hufken~Date,type='l',lwd=3,col=4))

# 
with(hufken.pace.pred.hghs,
     points(Date[which.max(cover.hufken)],cover.hufken[which.max(cover.hufken)],cex=2,pch=3))
with(hufken.pace.pred.hgls,
     points(Date[which.max(cover.hufken)],cover.hufken[which.max(cover.hufken)],cex=2,pch=3))
with(hufken.pace.pred.lgls,
     points(Date[which.max(cover.hufken)],cover.hufken[which.max(cover.hufken)],cex=2,pch=3))
with(hufken.pace.pred.lghs,
     points(Date[which.max(cover.hufken)],cover.hufken[which.max(cover.hufken)],cex=2,pch=3))

# 
legend('topright',legend = c('Onset of Brown-down',
                             'Low growth sensitivity; Low senescence sensitivity', 
                             'Low growth sensitivity; High senescence sensitivity',
                             'High growth sensitivity; Low senescence sensitivity',
                             'High growth sensitivity; High senescence sensitivity'),bty='n',
       pch=c(3,rep(NA,4)),
       lty = 1,lwd=3,
       col = c('grey',palette()))
abline(v=hufken.pace.pred.lgls$Date[20],lty='dashed',col='navy')
legend('topleft',legend = '(a)',bty='n')
# b#########
# first try one set of parameter
swc.capacity <- 0.3
swc.wilt <- 0.05
cover.max <- 1
cover.in <- NA
cover.vec <- seq(0,1,by=0.01)
#ll
q <- 0.6
q.s <- 4

# calculate threshold
threshold.vec.ll<- c()

for (i in seq_along(cover.vec)){
  cover.in <- cover.vec[i]
  
  threshold.vec.ll[i] <- uniroot(solve.intersect.func,interval = c(0.05,0.3))$root
}
#
q <- 0.6
q.s <- 0.1

# calculate threshold
threshold.vec.lh<- c()

for (i in seq_along(cover.vec)){
  cover.in <- cover.vec[i]
  
  threshold.vec.lh[i] <- uniroot(solve.intersect.func,interval = c(0.05,0.3))$root
}
#
q <- 1
q.s <- 4

# calculate threshold
threshold.vec.hl<- c()

for (i in seq_along(cover.vec)){
  cover.in <- cover.vec[i]
  
  threshold.vec.hl[i] <- uniroot(solve.intersect.func,interval = c(0.05,0.3))$root
}
#
q <- 1
q.s <- 0.1

# calculate threshold
threshold.vec.hh<- c()

for (i in seq_along(cover.vec)){
  cover.in <- cover.vec[i]
  
  threshold.vec.hh[i] <- uniroot(solve.intersect.func,interval = c(0.05,0.3))$root
}

# # try another set of sensitivity
# q <- 0.5
# q.s <- 2
# cover.vec <- seq(0,1,by=0.01)
# threshold.vec.tmp.2 <- c()
# 
# for (i in seq_along(cover.vec)){
#   cover.in <- cover.vec[i]
#   
#   threshold.vec.tmp.2[i] <- uniroot(solve.intersect.func,interval = c(0.05,0.3))$root
# }
# # make plot
# # load color

# 
plot(threshold.vec.ll~cover.vec,ylim=c(0.05,0.3),
     xlab = 'Cover',ylab = 'Soil moisture',
     type='l',lwd=2,col = 1)
points(threshold.vec.lh~cover.vec,
       type='l',lwd=2,col = 2)
points(threshold.vec.hl~cover.vec,
       type='l',lwd=2,col = 3)
points(threshold.vec.hh~cover.vec,
       type='l',lwd=2,col = 4)
legend('topleft',legend = '(b)',bty='n')

dev.off()