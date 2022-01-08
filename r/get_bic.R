# #######################################
# get stats for fitting##
# #######################################

# function to get BIC value for each model#####
get.bic.func <- function(model.vec,data.vec,n.fd){
  # old way of doing bic based on assumed ll########
  # d.vec <- model.vec
  # m.vec <- data.vec
  # 
  # res <- model.vec - data.vec
  # n <- length(data.vec)  
  # w <- rep(1,n) #not applicable
  # 
  # ll <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  # # Cll-logLik(m)==0 #TRUE
  # k.original<-n.fd
  # df.ll<-k.original+1 
  # bic<- -2 * ll + log(n) * df.ll
  
  # new way to do bic bsaed on mse####
  # https://machinelearningmastery.com/probabilistic-model-selection-measures/
  n <- length(data.vec)  
  mse <- mean((model.vec - data.vec)^2,na.rm=T)
  # BIC = -2 * LL + log(N) * k
  bic <- n * log(mse) + n.fd * log(n)
  
  return(bic)
}
# 
rmse.func = function(m, o){
  sqrt(mean((m - o)^2,na.rm=T))
}

# 
get.r.func <- function(x.vec,y.vec){
  fit.lm <- lm(y.vec~x.vec)
  return(summary(fit.lm)$r.squared)
}

# get fit####
source('r/read_spc_nm.R')

# process ym data and pace data
pace.ls <- list()

for (i in seq_along(species.vec)) {
  fn.con.11 <- sprintf('tmp/pred.smsmv13.2q.chain.%s.Control.Ambient.rds',
                    species.vec[i])
  
  dat.con.11 <- readRDS(fn.con.11)
  dat.con.11 <- dat.con.11[,c('cover','cover.hufken')]

  pace.ls[[i]] <-  dat.con.11
  
}

pace.effect.df <- do.call(rbind,pace.ls)
pace.effect.df <- pace.effect.df[complete.cases(pace.effect.df),]

# same thing but for v10 process
pace.ls.v10 <- list()

for (i in seq_along(species.vec)) {
  fn.con.10 <- sprintf('tmp/pred.smv13.q1.qs0.chain.%s.Control.Ambient.rds',
                    species.vec[i])
  
  dat.con.10 <- readRDS(fn.con.10)
  dat.con.10 <- dat.con.10[,c('cover','cover.hufken')]
  
  pace.ls.v10[[i]] <-  dat.con.10
  
}

pace.effect.df.v10 <- do.call(rbind,pace.ls.v10)
pace.effect.df.v10 <- pace.effect.df.v10[complete.cases(pace.effect.df.v10),]

# get bic####
# v11
get.bic.func(pace.effect.df$cover.hufken,
             pace.effect.df$cover,n.fd=5)
get.r.func(pace.effect.df$cover.hufken,
           pace.effect.df$cover)
rmse.func(pace.effect.df$cover.hufken,
          pace.effect.df$cover)
# v10
get.bic.func(pace.effect.df.v10$cover.hufken,
             pace.effect.df.v10$cover,n.fd=3)
get.r.func(pace.effect.df.v10$cover.hufken,
           pace.effect.df.v10$cover)
rmse.func(pace.effect.df.v10$cover.hufken,
          pace.effect.df.v10$cover)

# plot(pace.effect.df.v10$cover.hufken~pace.effect.df$cover.hufken)
# pace.effect.df.v10$cover.hufken-pace.effect.df$cover.hufken
