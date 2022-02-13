source('r/read_spc_nm.R')
source('r/load.R')

out.df <- data.frame(spc = species.vec,
                     f.t.opt =NA,
                     f.extract =NA,
                     f.sec =NA,
                     f.growth = NA,
                     q = NA,
                     q.s = NA)
for(i in seq_along(species.vec)){
  
  fn <- sprintf('cache/v13.2q.chain.%s.bestfit.rds',species.vec[i])
  print(paste0('par file used: ',fn))

  chain.fes =  readRDS(fn)
  # 
  fit.par.vec <- subset(chain.fes[which(chain.fes$ll==max(chain.fes$ll)),],select=-c(ll))
  # 
  out.df$q[i] <- fit.par.vec$q
  out.df$q.s[i] <- fit.par.vec$q.s

}


beta.func <- function(x,a=0.05,b=0.3,q=5,is.q.s = FALSE){
 if(is.q.s){
   (1-x)^q
 }else{
   (x)^q
 }
}

swc.vec <- seq(0,1,by=0.001)
beta.growth.ls <- beta.sene.ls <- list()
for (i.nm in seq_along(species.vec)) {
  beta.growth.ls[[i.nm]] <- beta.func(x=swc.vec,q=out.df$q[i.nm])
  beta.sene.ls[[i.nm]] <- beta.func(x=swc.vec,q=out.df$q.s[i.nm],is.q.s = TRUE)
}



png('figures/fw.png',width = 600,height = 2*600*0.618)
par(mar=c(5,5,1,1),mfrow=c(2,1))
palette(c(col.df$iris))

col.nm.vec <- c(1,1,2,2,2,3,3,3,4,4)
lty.vec <- c(1,2,1,2,3,1,2,3,1,2)
plot(beta.growth.ls[[1]]~swc.vec,type='l',col=col.nm.vec[1],
     xlab='Soil moisture availability',ylab=expression(beta[growth]),lwd=3,lty=lty.vec[1])
# 
for (i in seq_along(beta.growth.ls)) {
  points(beta.growth.ls[[i]]~swc.vec,type='l',col=col.nm.vec[i],lty=lty.vec[i],lwd=3)
}
legend('topleft',legend = '(a)',bty='n')


legend('bottomright',legend = species.vec.nm,
       lty=lty.vec,col=col.nm.vec,ncol = 2,lwd=3)

# 
plot(beta.sene.ls[[1]]~swc.vec,type='l',col=col.nm.vec[1],xlim=c(0,1),
      xlab='Soil moisture availability',ylab=expression(beta[senescence]),lty=lty.vec[1],lwd=3)


for (i in seq_along(beta.sene.ls)) {
  points(beta.sene.ls[[i]]~swc.vec,type='l',col=col.nm.vec[i],lty=lty.vec[i],lwd=3)
}
legend('topleft',legend = '(b)',bty='n')
dev.off()

