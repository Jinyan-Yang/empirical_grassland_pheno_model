
t.mean <- seq(17,45,by=0.1)

# f.t.opt <- 15

t.func <- function(t.mean,f.t.opt,t.max){
  # need to determine which is the best way#
  
  # #1.ch tyope of t dependence
  # return((t.max-t.mean)/(t.max-f.t.opt)*(t.mean/f.t.opt)^(f.t.opt/(t.max-f.t.opt)))
  
  # #2. p based sharp decline
  # h.val <- pnorm(t.mean,mean=f.t.opt,sd = 5,lower.tail = F)
  # l.val <- pnorm(t.mean,mean=f.t.opt,sd = 5,lower.tail = T)
  # return(min(h.val,l.val)*2)
  
  
  # 3. sgs
  # frac <- ((t.mean - t.mn)/(t.r - t.mn))^q *
  #   ((1+q)*f.t.opt - t.mn - q*t.mean) / 
  #   ((1+q)*f.t.opt - t.mn - q*t.r)

  q <- 2.5#power or sensitivity of the function
  t.mn <- 5#MIN t FOR GROWTH; 5 degree seems reasonable?
  
  # here we assume that Topt and 
  # Tr (the T at which function is 1) to be the same 
  # so that the function is from 0 to 1
  frac <- ((t.mean - t.mn)/(f.t.opt - t.mn))^q *
    ((1+q)*f.t.opt - t.mn - q*t.mean) /
    ((1+q)*f.t.opt - t.mn - q*f.t.opt)
  
  return(frac)

}

# #####testing T func####
# f.vec.35 <- c()
# f.vec.20 <- c()
# for(i in seq_along(t.mean)){
#   f.vec.35[i] <- t.func(t.mean = t.mean[i],
#                      f.t.opt = 35,t.max=45)
#   f.vec.20[i] <- t.func(t.mean = t.mean[i],
#                         f.t.opt = 20,t.max=45)
# }

# plot(f.vec.20~t.mean,ylim=c(0,1),col='navy')
# points(f.vec.35~t.mean,col='red')

