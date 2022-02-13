# ###############################
# calculate onset of browndown
# ###############################

out.sample.ls <- list()

# Loop over species: randomly sample parameter values from the output chains
for(i in  seq_along(species.vec)){
  fn.1 <- sprintf('cache/v13.2q.chain.%s.bestfit.rds',species.vec[i])
  chain.fes <- readRDS(fn.1)
  chain.fes$spc <- species.vec[i]
  # make a list of all species 
  out.sample.ls[[i]] <- chain.fes
}

# This function finds the dV/dt for a given soil moisture content
# Not actually needed! 
# Also, ideally would be passing parameters to this
solve.intersect.func <- function(swc.norm){
  cover.max <- 1
  loss.f <- swc.norm^q
  loss.f.s <- (1-swc.norm)^q.s
  growth.vec <- f.growth* loss.f * (1 - cover.in / cover.max)
  senescence.vec <- f.sen*(loss.f.s) * cover.in
  return(growth.vec - senescence.vec)
}

# This function finds the threshold value of Vt for which dV/dt = 0 at a given SWC
vthresh <- function(swc, rg, qg, rs, qs) {
  rg*swc^qg / (rs*(1-swc)^qs + rg*swc^qg)
}


# Make a list to hold CI's
pred.ci.ls <- list()

# species.vec is the list of species
for (iter.nm in seq_along(species.vec)) {
  
  # This is the 200 random samples of the 6 parameters for this species
  tmp.df <- out.sample.ls[[iter.nm]]
  
  # Letting SWC vary from 0 to 1 in 0.01 steps
  swc.vec <- seq(0,1,by=0.005)
  
  # calculate threshold
  
  # Set up df to hold values. Length 101 (corresponds to swc values) and columns are quantiles 
  tmp.x.df <- data.frame(swc=rep(NA,length(swc.vec)),
                         q.05=NA,
                         q.5=NA,
                         q.95=NA,
                         spc = species.vec[iter.nm])
  
  # For each value of swc
  for (i in seq_along(swc.vec)){
    threshold.vec.tmp <- c()
    swc.in <- swc.vec[i]
    
    # For each random sample of parameters
    for(j in seq_along(tmp.df$spc)){
      f.sen <- tmp.df[j,3]
      f.growth <- tmp.df[j,4]
      q <- tmp.df[j,5]
      q.s <- tmp.df[j,6]
      
      # Find the soil moisture value for which dV/dt = 0
      x <- vthresh(swc = swc.in, rg = f.growth, qg = q, rs = f.sen, qs = q.s)
      
      # append the solution to the end of the vector
      threshold.vec.tmp <- c(threshold.vec.tmp, x)
    }
    
    # Find the quantiles of the soil moisture values given the variation in parameters
    q.tmp <-  quantile(threshold.vec.tmp,probs = c(.05,.5,.95),na.rm=T)
    # and add quantiles to the data frame
    tmp.x.df$swc[i] <-  swc.in
    tmp.x.df$q.05[i] <- q.tmp[[1]]
    tmp.x.df$q.5[i] <- q.tmp[[2]]
    tmp.x.df$q.95[i] <- q.tmp[[3]]
  }
  # add quantiles to data frame
  pred.ci.ls[[iter.nm]] <- tmp.x.df
}

ci.df <- do.call(rbind,pred.ci.ls)
# saveRDS(ci.df,'tmp/threshold_BM.rds')
# 
palette(c(col.df$iris))


png('figures/threshold.png',width = 800,height = 800)
layout(matrix(c(1:6,7,8,11,9,10,11),3,4, byrow = FALSE), 
       # c(3,1), c(1,3),
       respect = TRUE)

par(mar=rep(2,4),xpd=TRUE,oma=rep(2.2,4))



# par(mfrow=c(4,3))
letter.nm=1
for (i in c(3:10,1,2)) {
  plot(0,pch=16,col='white',xlim=c(0,1),ylim=c(0,1),
       xlab=' ',ylab=' ')
  
  mid.df <- pred.ci.ls[[i]]$q.5
  
  polygon(c(0,swc.vec,1),c(0,mid.df,0),col=col.df$auLandscape[2])
  polygon(c(0,1,rev(swc.vec),0),c(1,1,rev(mid.df),0),col=col.df$auLandscape[3])
  # species.vec.nm <- species.vec
  # species.vec.nm[species.vec.nm=='Flux'] <- 'Flux Tower'
  legend('topleft',legend = paste0('(',letters[letter.nm],') ',species.vec.nm[i]))
  
  points(q.05~swc,data = pred.ci.ls[[i]],type='l',lty='dashed',lwd=2)
  points(q.95~swc,data = pred.ci.ls[[i]],type='l',lty='dashed',lwd=2)
  if(letter.nm==1){
    legend('bottomright',
           legend = c('Green-up','Brown-down'),
           pch=15,col=col.df$auLandscape[2:3])
  }
  if(letter.nm==2){
    mtext('Cover',side=2,xpd=T,line=3)
  }
  if(letter.nm==6){
    mtext('Soil moisture availability',side=1,xpd=T,line=3,adj=1)
  }
  letter.nm =letter.nm+1
}
# par(mar=c(5,5,3,1))

plot(q.5~swc,
     data = pred.ci.ls[[10]],
     ylim=c(0,1),
     xlab = ' ',
     ylab = ' ',
     type='l',lwd=2,col = 4,lty=2)

points(q.5~swc,data = pred.ci.ls[[1]],type='l',col=1,lwd=2,lty=1)
points(q.5~swc,data = pred.ci.ls[[2]],type='l',col=1,lwd=2,lty=2)
points(q.5~swc,data = pred.ci.ls[[3]],type='l',col=2,lwd=2,lty=1)
points(q.5~swc,data = pred.ci.ls[[4]],type='l',col=2,lwd=2,lty=2)
points(q.5~swc,data = pred.ci.ls[[5]],type='l',col=2,lwd=2,lty=3)
points(q.5~swc,data = pred.ci.ls[[6]],type='l',col=3,lwd=2,lty=1)
points(q.5~swc,data = pred.ci.ls[[7]],type='l',col=3,lwd=2,lty=2)
points(q.5~swc,data = pred.ci.ls[[8]],type='l',col=3,lwd=2,lty=3)
points(q.5~swc,data = pred.ci.ls[[9]],type='l',col=4,lwd=2,lty=1)

legend('topleft',legend = '(k)',bty='n')
# species.vec.nm <- species.vec
# species.vec.nm[species.vec.nm=='Flux'] <-' Flux Tower'
legend('bottomright',legend = species.vec.nm,
       lty=c(1,2,1,2,3,1,2,3,1,2),col=c(1,1,2,2,2,3,3,3,4,4),lwd=2,
       ncol = 2)



dev.off()

# matrix(c(1:6,7,8,11,9,10,11),3,4, byrow = FALSE)
# matrix(c(3:10,11,1,2,11),3,4, byrow = FALSE)

