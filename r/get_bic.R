# function to get BIC value for each model#####
get.bic.func <- function(model.vec,data.vec,n.fd){
  d.vec <- model.vec
  m.vec <- data.vec
  
  res <- model.vec - data.vec
  n <- length(data.vec)  
  w <- rep(1,n) #not applicable
  
  ll<-0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  # Cll-logLik(m)==0 #TRUE
  k.original<-n.fd
  df.ll<-k.original+1 
  bic<- -2 * ll + log(n) * df.ll
  
  return(bic)
  # -0.5 *  sum((m.vec - d.vec)^2)/d.sd - n.fd*log(length(m.vec))
}


get.r.func <- function(x.vec,y.vec){
  # x.vec = optbb.spots$ALEAF
  # y.vec=spot.amb$Photo
  fit.lm <- lm(y.vec~x.vec)
  return(summary(fit.lm)$r.squared)
}

# get fit####
source('r/read_spc_nm.R')

# process ym data and pace data
pace.ls <- list()

for (i in seq_along(species.vec)) {
  fn.con <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.Ambient.rds',
                    species.vec[i])
  
  dat.con <- readRDS(fn.con)
  dat.con <- dat.con[,c('cover','cover.hufken')]

  pace.ls[[i]] <-  dat.con
  
}

pace.effect.df <- do.call(rbind,pace.ls)
pace.effect.df <- pace.effect.df[complete.cases(pace.effect.df),]

# get bic####
get.bic.func(pace.effect.df$cover.hufken,pace.effect.df$cover,n.fd=5)
get.r.func(pace.effect.df$cover.hufken,pace.effect.df$cover)
