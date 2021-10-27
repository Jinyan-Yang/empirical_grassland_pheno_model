
# species.vec <- c('Luc','Fes','Rye','Kan',
#                  'YM')
# species.vec <- c("Bis",    "Dig",  "Fes",    "Kan",    
#                  "Luc",  "Rho",    "Rye",
#                  'YM','Flux')
species.vec <- c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM','Flux')

out.df <- data.frame(spc = species.vec,
                     f.t.opt =NA,
                     f.extract =NA,
                     f.sec =NA,
                     f.growth = NA,
                     q = NA,
                     q.s = NA,
                     q.05=NA,
                     q.95=NA,
                     q.s.05=NA,
                     q.s.95=NA)
for(i in seq_along(species.vec)){
  fn <- sprintf('cache/smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  
  chain.3.ls = readRDS(fn)
  
  chain.3.ls.new = lapply(chain.3.ls,function(m.in)m.in[round(2*nrow(m.in)/3):nrow(m.in),])
  chain.fes <- do.call(rbind,chain.3.ls.new)
  
  
  fitted.val <- colMeans(chain.fes)
  
  out.df[i,2:7] <- fitted.val
  
  q.quant <- quantile(chain.fes[,5],probs = c(.05,.95))
  qs.quant <- quantile(chain.fes[,6],probs = c(.05,.95))
  out.df$q.05[i] <- q.quant[[1]]
  out.df$q.95[i] <- q.quant[[2]]
  
  out.df$q.s.05[i] <- qs.quant[[1]]
  out.df$q.s.95[i] <- qs.quant[[2]]
}


beta.func <- function(x,a=0.05,b=0.3,q=5,is.q.s = FALSE){
 if(is.q.s){
   (1-(x - a) / (b - a))^q 
 }else{
   ((x - a) / (b - a))^q 
 }
}

swc.vec <- seq(0.05,0.3,by=0.001)
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
     xlab='Soil moisture',ylab=expression(f[w]),lwd=3,lty=lty.vec[1])

# points(beta.growth.ls[[2]]~swc.vec,type='l',col=2,lwd=3)
# points(beta.growth.ls[[3]]~swc.vec,type='l',col=3,lwd=3)
# points(beta.growth.ls[[4]]~swc.vec,type='l',col=4,lwd=3)
# points(beta.growth.ls[[5]]~swc.vec,type='l',col=5,lwd=3)


for (i in seq_along(beta.growth.ls)) {
  points(beta.growth.ls[[i]]~swc.vec,type='l',col=col.nm.vec[i],lty=lty.vec[i],lwd=3)
}
legend('topleft',legend = '(a)',bty='n')
out.df$spc[out.df$spc=='Flux'] <- 'Flux Tower'
legend('bottomright',legend = out.df$spc,
       lty=lty.vec,col=col.nm.vec,ncol = 2,lwd=3)

# 
plot(beta.sene.ls[[1]]~swc.vec,type='l',col=col.nm.vec[1],xlim=c(0.03,0.3),
      xlab='Soil moisture',ylab=expression(f[w*','*s]),lty=lty.vec[1],lwd=3)

# points(beta.sene.ls[[2]]~swc.vec,type='l',col=2,lwd=3,lty='dashed')
for (i in seq_along(beta.sene.ls)) {
  points(beta.sene.ls[[i]]~swc.vec,type='l',col=col.nm.vec[i],lty=lty.vec[i],lwd=3)
}
legend('topleft',legend = '(b)',bty='n')
dev.off()

