source('r/plot.mcmc.r')
source('r/read_spc_nm.R')
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
# # make PDF plots#####
# pdf('figures/v13.2q.pace.pdf',width = 8,height = 8*0.618)
# 
# # ym
# plot.mcmc.func.2q(ym.18.df,species.in='ym',prep.in = 'Control','Ambient',
#                   my.fun = phenoGrass.func.v13,
#                   nm.note='v13.2q',use.smooth = TRUE,swc.in.cap = 0.3,swc.in.wilt = 0.05,bucket.size = 1000)
# plot.title.func('YM') 
# # pace
# for (i in c(1,2,3,4,6,10)) {
#   plot.mcmc.func.2q(gcc.met.pace.df,species.vec[i],'Control','Ambient',
#                     my.fun = phenoGrass.func.v13,
#                     nm.note='v13.2q',use.smooth = TRUE,day.lag = 3)
#   plot.title.func(species.vec[i]) 
#   
# }
# # # cw
# # for(i in seq_along(site.vec[1:7])){
# #   fn <- sprintf('cache/smv13chain.%s.Control.Ambient.rds',site.vec[i])
# #   if(file.exists(fn)){
# # 
# #     plot.mcmc.func.2q(gcc.met.cw.df,site.vec[i],'Control','Ambient'
# #                       ,subplot = NULL,nm.note = 'v13.2q',use.smooth = TRUE,my.fun =phenoGrass.func.v13,
# #                    swc.in.wilt = 0.05,swc.in.cap = 0.3,bucket.size=1000)
# #     plot.title.func(site.vec[i])
# #   }
# # }
# 
# dev.off()


# pdf('figures/v13.2q.diag.pdf',width = 8,height = 8*0.618)
# fit pace####

# species.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','YM','Flux')
# species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM','Flux')
# species.vec.nm <- c('Bis','Med','Dig','The','Chl','Fes','Pha','Lol','YM','Flux Tower')
# species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM','Flux')
# legume = 1; C4=2,C3=3,mix=4
# c(1,1,2,2,2,3,3,3,4,4)
# out.df <- data.frame(spc = species.vec,
#                      f.t.opt =NA,
#                      f.extract =NA,
#                      f.sec =NA,
#                      f.growth = NA,
#                      q = NA,
#                      q.s = NA)
# for(i in  seq_along(species.vec)){
#   fn <- sprintf('cache/smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
#   
#   chain.3.ls = readRDS(fn)
#   
#   chain.3.ls.new = lapply(chain.3.ls,function(m.in)m.in[round(2*nrow(m.in)/3):nrow(m.in),])
#   chain.fes <- do.call(rbind,chain.3.ls.new)
#   
#  
#   fitted.val <- colMeans(chain.fes[40000:nrow(chain.fes),])
#   
#   out.df[i,2:7] <- fitted.val
#   
#   # lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i])
#   
#   # par(mfrow=c(3,2),mar=c(5,5,1,1))
#   # for(par.num in 1:6){
#   #   
#   #   plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(15000:nrow(chain.3.ls[[1]])))
#   #   
#   # }
# }

out.sample.ls <- list()

for(i in  seq_along(species.vec)){
  fn <- sprintf('cache/smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  
  chain.3.ls = readRDS(fn)
  
  chain.3.ls.new = lapply(chain.3.ls,function(m.in)m.in[round(nrow(m.in)*0.75):nrow(m.in),])
  chain.fes <- do.call(rbind,chain.3.ls.new)
  
  # 
  set.seed(1935)
  sample.index <- sample(1:nrow(chain.fes),200)
  
  # tmp.m <- apply(chain.fes, 2, sample,size=100)
  
  tmp.m <- as.data.frame(chain.fes[sample.index,])
  tmp.m <- as.data.frame(tmp.m)
  tmp.m$spc <- species.vec[i]
  
  
  out.sample.ls[[i]] <- tmp.m
  
  # lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i])
  
  # par(mfrow=c(3,2),mar=c(5,5,1,1))
  # for(par.num in 1:6){
  #   
  #   plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(15000:nrow(chain.3.ls[[1]])))
  #   
  # }
}
# hist(tmp.m[,2])


# get.threshold.func <- function(gcc.vec = seq(0.1,1,by=0.01),
#                                f.sec =0.25,
#                                f.growth = 0.095,
#                                swc.wilt =0.05 ,
#                                swc.capacity = 0.3,
#                                bucket.size=1000,
#                                q = 0.05,
#                                q.s = 1.5){
#   # gcc.vec <- seq(0.1,1,by=0.01)
#   # swc.vec <- seq(0.05,0.3,by=0.01)
#   pred.ls <- list()
#   break.points <- c()
#   
#   for (i in seq_along(gcc.vec)) {
#     
#     # set rain scenario
#     tmp.df$vwc <- 0.3
#     tmp.df$GCC.norm.smooth <- gcc.vec[i]
#     
#     pred.out <- phenoGrass.func.v13(tmp.df,
#                                     f.h = 200,
#                                     f.t.opt = 30,
#                                     f.extract = 0.0287,
#                                     f.sec =f.sec,
#                                     f.growth = f.growth,
#                                     swc.wilt =swc.wilt ,
#                                     swc.capacity = swc.capacity,
#                                     bucket.size=bucket.size,
#                                     t.max = 45,
#                                     day.lay = 1,
#                                     use.smooth=TRUE,
#                                     q = q,
#                                     q.s = q.s)
#     pred.out <- pred.out[complete.cases(pred.out),]
#     pred.out$net.change <- pred.out$growth - pred.out$senescence
#     
#     # library(mgcv)
#     # 
#     # fit.gam <- gam(vwc.hufken~s(net.change,k=4),data = pred.out)
#     # 
#     # predict(fit.gam,newdata = data.frame(net.change=0))
#     
#     break.points[i] <- pred.out$vwc.hufken[abs(pred.out$net.change)<0.05]
#     
#     # break.points[i] <- pred.out$vwc.hufken[pred.out$senescence ==
#     #                                          max(pred.out$senescence,
#     #                                              na.rm=T)]
#   }
#   
#   threshold.df <- data.frame(gcc = gcc.vec,
#                              threshold = break.points) 
#   
#   return(threshold.df)
# }
# 
# out.df$swc.capacity <- 0.3#c(0.13,0.13,0.13,0.3)
# out.df$bucket.size  <- 1000#c(300,300,300,1000)
# 
# pred.ls <- list()
# 
# for (iter.nm in 1:nrow(out.df)) {
#   pred.ls[[iter.nm]] <- get.threshold.func(f.sec=out.df$f.sec[iter.nm],
#                                            f.growth = out.df$f.growth[iter.nm],
#                                            q = out.df$q[iter.nm],
#                                            q.s = out.df$q.s[iter.nm],
#                                            swc.wilt =0.05 ,
#                                            swc.capacity = out.df$swc.capacity[iter.nm],
#                                            bucket.size=out.df$bucket.size[iter.nm])
#   
#   pred.ls[[iter.nm]]$spc <- out.df$spc[iter.nm]
# }
# 
# 
# palette(col.df$auLandscape)
# # same plot with fitted values from ym and pace
# plot(threshold~gcc,data = pred.ls[[1]],type='l',
#      xlab='Initial cover',
#      ylab='SWC threshold',ylim=c(0,0.3),col=1)
# points(threshold~gcc,data = pred.ls[[2]],type='l',col=2)
# points(threshold~gcc,data = pred.ls[[3]],type='l',col=3)
# points(threshold~gcc,data = pred.ls[[4]],type='l',col=4)
# abline(h=0.05,lty='dashed')
# legend('topleft',legend = out.df$spc,
#        lty='solid',col=palette())
solve.intersect.func <- function(swc.norm){
  # swc.norm <-  (vwc.in- swc.wilt) / (swc.capacity - swc.wilt)
  loss.f <- swc.norm^q
  loss.f.s <- (1-swc.norm)^q.s
  
  growth.vec <- f.growth*
    loss.f *
    (1 - cover.in / cover.max)
  
  senescence.vec <- f.sen*(loss.f.s) *
    # (1 - cover.pred.vec[nm.day-1])*
    cover.in
  
  return(growth.vec - senescence.vec)
}

# find.thresh.func <- function(q,q.s){
#   swc.capacity <- 0.3
#   swc.wilt <- 0.05
# 
#   cover.max <- 1
#   cover.in <- 0.5
#   
#   cover.vec <- seq(0,1,by=0.01)
#   # calculate threshold
#   threshold.vec.tmp <- c()
#   
#   for (i in seq_along(cover.vec)){
#     cover.in <- cover.vec[i]
#     
#     threshold.vec.tmp[i] <- uniroot(solve.intersect.func,interval = c(0.05,0.3),q =q,q.s=q.s)$root
#     
#   }
#   # print(threshold.vec.tmp)
#   return(threshold.vec.tmp)
# }
# 
# # find.thresh.func(q=2,q.s=0.5)
# 
# pred.ls <- list()
# for (iter.nm in 1:nrow(out.df)) {
#   # 
#   if(species.vec[i]=='ym'){
# 
#     ym.met.df <- readRDS('cache/ym/ym.met.rds')
#     # 
#     swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
#     swc.cap = round(swc.ym.con[[2]]*10)/10
#     swc.wilt = round(swc.ym.con[[1]]*100)/100
# 
#   }else if(species.vec[i]=='flux'){
# 
# 
#     swc.cap =  0.4
#     swc.wilt = 0.01
# 
#   }else{
#     swc.cap = 0.13
#     swc.wilt = 0.05
#     bucket.size=300
#   }
# 
#   q <- out.df$q[iter.nm]
#   q.s <- out.df$q.s[iter.nm]
#   f.growth <- out.df$f.growth[iter.nm]
#   f.sen <- out.df$f.sec[iter.nm]
#   
#   cover.max <- 1
#   cover.in <- 0.5
#   
#   cover.vec <- seq(0,1,by=0.01)
#   
#   # calculate threshold
#   threshold.vec.tmp <- c()
#   
#   for (i in seq_along(cover.vec)){
#     cover.in <- cover.vec[i]
#     
#     threshold.vec.tmp[i] <- uniroot(solve.intersect.func,interval = c(swc.wilt,swc.cap))$root
#   }
#   
#   
#   pred.ls[[iter.nm]] <- threshold.vec.tmp
# }
# get ci
pred.ci.ls <- list()

for (iter.nm in seq_along(species.vec)) {
    # if(species.vec[iter.nm]=='YM'){
    # 
    #   ym.met.df <- readRDS('cache/ym/ym.met.rds')
    #   #
    #   swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    #   swc.cap = round(swc.ym.con[[2]]*10)/10
    #   swc.wilt = round(swc.ym.con[[1]]*100)/100
    # 
    # }else if(species.vec[i]=='Flux'){
    #   swc.cap =  0.4
    #   swc.wilt = 0.01
    # 
    # }else{
    #   swc.cap = 0.13
    #   swc.wilt = 0.05
    #   bucket.size=300
    # }
  # # 
  # swc.capacity <- 0
  # swc.wilt <- 1
  # 
  tmp.df <- out.sample.ls[[iter.nm]]

  cover.max <- 1
  # cover.in <- 0.5
  
  cover.vec <- seq(0,1,by=0.01)
  
  # calculate threshold
 
  tmp.x.df <- data.frame(cover=rep(NA,length(cover.vec)),
                         q.05=NA,
                         q.5=NA,
                         q.95=NA,
                         spc = species.vec[iter.nm])
  for (i in seq_along(cover.vec)){
    threshold.vec.tmp <- c()
    cover.in <- cover.vec[i]
    for(j in seq_along(tmp.df$spc)){
      q <- tmp.df[j,5]
      q.s <- tmp.df[j,6]
      f.growth <- tmp.df[j,3]
      f.sen <- tmp.df[j,4]
      
      x <- try(uniroot(solve.intersect.func,interval = c(0,1))$root)
      
      if(class(x)=='try-error'){
        x = NA
      }
      
      threshold.vec.tmp <- c(threshold.vec.tmp, x)
    }
    q.tmp <-  quantile(threshold.vec.tmp,probs = c(.05,.5,.95),na.rm=T)
    
    tmp.x.df$cover[i] <-  cover.in
    tmp.x.df$q.05[i] <- q.tmp[[1]]
    tmp.x.df$q.5[i] <- q.tmp[[2]]
    tmp.x.df$q.95[i] <- q.tmp[[3]]
}
  
  pred.ci.ls[[iter.nm]] <- tmp.x.df
}
ci.df <- do.call(rbind,pred.ci.ls)
# 
palette(c(col.df$iris))


# # species.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','ym')
# plot(0,pch=16,col='white',xlim=c(0,1),ylim=c(0,0.3),
#      xlab='Cover',ylab='Soil moisture')
# polygon(c(0,cover.vec,1),c(0,pred.ls[[1]],0),col=col.df$auLandscape[3])
# polygon(c(0,1,rev(cover.vec),0),c(0.3,0.3,rev(pred.ls[[1]]),0),col=col.df$auLandscape[2])
# legend(x=0.85,y=0.05,legend = c('Green-up','Brown-down'),pch=15,col=col.df$auLandscape[2:3],)


# png('figures/threshold.png',width = 600,height = 600*0.618)
# plot(pred.ls[[10]]~cover.vec,
#      ylim=c(0.05,0.3),
#      xlab = 'Cover',
#      ylab = 'Soil moisture',
#      type='l',lwd=2,col = 4,lty=1)
# 
# points(pred.ls[[1]]~cover.vec,type='l',col=1,lwd=2,lty=1)
# points(pred.ls[[2]]~cover.vec,type='l',col=1,lwd=2,lty=2)
# points(pred.ls[[3]]~cover.vec,type='l',col=2,lwd=2,lty=1)
# points(pred.ls[[4]]~cover.vec,type='l',col=2,lwd=2,lty=2)
# points(pred.ls[[5]]~cover.vec,type='l',col=2,lwd=2,lty=3)
# points(pred.ls[[6]]~cover.vec,type='l',col=3,lwd=2,lty=1)
# points(pred.ls[[7]]~cover.vec,type='l',col=3,lwd=2,lty=2)
# points(pred.ls[[8]]~cover.vec,type='l',col=3,lwd=2,lty=3)
# 
# points(pred.ls[[9]]~cover.vec,type='l',col=4,lwd=2,lty=2)
# 
# legend('topleft',legend = '(a)',bty='n')
# legend('bottomright',legend = out.df$spc,
#        lty=c(1,2,1,2,3,1,2,3,1,2),col=c(1,1,2,2,2,3,3,3,4,4),lwd=2,
#        ncol = 2)
# dev.off()

png('figures/threshold.png',width = 800,height = 800)
layout(matrix(c(1:6,7,8,11,9,10,11),3,4, byrow = FALSE), 
       # c(3,1), c(1,3),
       respect = TRUE)

par(mar=rep(2,4),xpd=TRUE,oma=rep(2,4))



# par(mfrow=c(4,3))
letter.nm=1
for (i in c(3:10,1,2)) {
  plot(0,pch=16,col='white',xlim=c(0,1),ylim=c(0,1),
       xlab=' ',ylab=' ')
  
  mid.df <- pred.ci.ls[[i]]$q.5
 
  polygon(c(0,cover.vec,1),c(0,mid.df,0),col=col.df$auLandscape[3])
  polygon(c(0,1,rev(cover.vec),0),c(1,1,rev(mid.df),0),col=col.df$auLandscape[2])
  # species.vec.nm <- species.vec
  # species.vec.nm[species.vec.nm=='Flux'] <- 'Flux Tower'
  legend('topleft',legend = paste0('(',letters[letter.nm],') ',species.vec.nm[i]))
  
  points(q.05~cover,data = pred.ci.ls[[i]],type='l',lty='dashed',lwd=2)
  points(q.95~cover,data = pred.ci.ls[[i]],type='l',lty='dashed',lwd=2)
  if(letter.nm==1){
    legend('bottomright',
           legend = c('Green-up','Brown-down'),
           pch=15,col=col.df$auLandscape[2:3])
  }
  if(letter.nm==2){
    mtext('Soil moisture availability',side=2,xpd=T,line=3)
  }
  if(letter.nm==6){
    mtext('Cover',side=1,xpd=T,line=3,adj=1)
  }
  letter.nm =letter.nm+1
}
# par(mar=c(5,5,3,1))

plot(q.5~cover,
     data = pred.ci.ls[[10]],
     ylim=c(0,1),
     xlab = ' ',
     ylab = ' ',
     type='l',lwd=2,col = 4,lty=2)

points(q.5~cover,data = pred.ci.ls[[1]],type='l',col=1,lwd=2,lty=1)
points(q.5~cover,data = pred.ci.ls[[2]],type='l',col=1,lwd=2,lty=2)
points(q.5~cover,data = pred.ci.ls[[3]],type='l',col=2,lwd=2,lty=1)
points(q.5~cover,data = pred.ci.ls[[4]],type='l',col=2,lwd=2,lty=2)
points(q.5~cover,data = pred.ci.ls[[5]],type='l',col=2,lwd=2,lty=3)
points(q.5~cover,data = pred.ci.ls[[6]],type='l',col=3,lwd=2,lty=1)
points(q.5~cover,data = pred.ci.ls[[7]],type='l',col=3,lwd=2,lty=2)
points(q.5~cover,data = pred.ci.ls[[8]],type='l',col=3,lwd=2,lty=3)

points(q.5~cover,data = pred.ci.ls[[9]],type='l',col=4,lwd=2,lty=1)

legend('topleft',legend = '(k)',bty='n')
# species.vec.nm <- species.vec
# species.vec.nm[species.vec.nm=='Flux'] <-' Flux Tower'
legend('bottomright',legend = species.vec.nm,
       lty=c(1,2,1,2,3,1,2,3,1,2),col=c(1,1,2,2,2,3,3,3,4,4),lwd=2,
       ncol = 2)



dev.off()

# matrix(c(1:6,7,8,11,9,10,11),3,4, byrow = FALSE)
# matrix(c(3:10,11,1,2,11),3,4, byrow = FALSE)

