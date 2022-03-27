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

pdf('figures/theory_demo.pdf',width = 10,height = 10*.62)
#####
par(mar=c(5,5,1,1))
par(mfrow=c(2,2))


# a####
beta.func <- function(x,q=5,is.q.s = FALSE){
  if(is.q.s){
    (1-x)^q 
  }else{
    (x)^q 
  }
}
# 
swc.vwc <- seq(0,1,by=0.01)
co.up.vec <- beta.func(swc.vwc,q=4)
co.down.vec <- beta.func(swc.vwc,q=0.25)

co.up.vec.s <- beta.func(swc.vwc,q=0.25,is.q.s=T)
co.down.vec.s <- beta.func(swc.vwc,q=4,is.q.s=T)

plot(co.up.vec~swc.vwc,type='l',lwd=3,col='red',
     xlab='Soil moisture availability',ylab='Fractional reduction of growth rate')
points(co.down.vec~swc.vwc,type='l',lwd=3,col='black')
abline(a=0,b=1,col='grey',lwd=3)
legend('topleft',legend = '(a)',bty='n')
# b
plot(co.up.vec.s~swc.vwc,type='l',lwd=3,col='red',lty='dashed',
     xlab='Soil moisture availability',ylab='Fractional reduction of senescence rate')
points(co.down.vec.s~swc.vwc,type='l',lwd=3,col='black',lty='dashed')

abline(a=1,b=-1,col='grey',lwd=3,lty='dashed')
legend('topleft',legend = '(b)',bty='n')
# c####
plot.shape.func <- function(x,a,b,c=3,d=15){
  return(
    ((b-x) / (b-a))^c * exp(c/d*(1-((b-x)/(b-a))^d))
  )
}

x.vec <- seq(0.01,22,by=0.1)
# 
y.vec.ll <- plot.shape.func(x.vec,
                            a = 2,b=10)
# which.min(abs(y.vec.ll - y.vec.hh))

plot(y.vec.ll~x.vec,ylim=c(0.1,1),type='l',lwd=2,col=1,
     xlim=c(0.7,8),
     xaxt='n',yaxt='n',xlab='Time',ylab = 'Cover')
# 
y.vec.hh <- plot.shape.func(x.vec,
                            a = 5,b=10,
                            c=3,d=2)

points(y.vec.hh~x.vec,ylim=c(0,1),type='l',lwd=2,col=2)

# 
points(x.vec[which.max(y.vec.ll)],y.vec.ll[which.max(y.vec.ll)],cex=2,pch=3,col='grey20')
points(x.vec[which.max(y.vec.hh)],y.vec.hh[which.max(y.vec.hh)],cex=2,pch=3,col='grey20')

# 
legend('bottom',legend = c('Onset of Brown-down'),bty='n',
       pch=c(3),
       col = c('grey20'))

legend('topleft',legend = '(c)',bty='n')

# #########
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


# 
plot(threshold.vec.ll~cover.vec,ylim=c(0.05,0.3),
     xlab = 'Cover',ylab = 'Soil moisture availability threshold',
     type='l',lwd=2,col = 'black',yaxt='n')
axis(2,at = seq(0,0.3,length.out = 5),labels = seq(0,1,length.out = 5))

points(threshold.vec.hh~cover.vec,
       type='l',lwd=2,col = 'red')
legend('topleft',legend = '(d)',bty='n')

dev.off()
