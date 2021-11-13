###################################################################
# functions to do the mcmc fitting fot V13 version of the model####
###################################################################

source('r/functions_mcmc_v12.r')
mh.MCMC.func <- function(iterations,par.df,
                         gcc.met.pace.df.16,
                         bucket.size,swc.wilt ,
                         swc.capacity,day.lay,
                         my.fun,
                         use.smooth){
  
  # get prior 
  prior.prob <- prior.func(par.df)
  
  # start MC ####
  # intial
  chain = array(dim = c(iterations+1,ncol(par.df)))
  chain[1,] = as.numeric(par.df['initial',])
  
  # chain move on
  for (i in 1:iterations){
    proposal = proposal.func(par.df)
    
    # prior.prob,data,data.sd,bucket.size = 300,...
    probab = exp(posterior.func(prior.prob,FUN = my.fun,
                                gcc.df = gcc.met.pace.df.16,
                                f.h = 222,
                                f.t.opt = proposal[1],
                                f.extract = proposal[2],
                                f.sec = proposal[3],
                                f.growth = proposal[4] ,
                                q = proposal[5] ,
                                t.max = 45,
                                day.lay = day.lay,
                                bucket.size = bucket.size,
                                swc.wilt = swc.wilt ,
                                swc.capacity = swc.capacity,
                                use.smooth = use.smooth) - 
                   posterior.func(prior.prob,FUN = my.fun,
                                  gcc.df = gcc.met.pace.df.16,
                                  f.h = 222,
                                  f.t.opt = chain[i,1],
                                  f.extract = chain[i,2],
                                  f.sec = chain[i,3],
                                  f.growth = chain[i,4] ,
                                  q = chain[i,5] ,
                                  t.max = 45,
                                  day.lay = day.lay,
                                  swc.wilt = swc.wilt ,
                                  swc.capacity = swc.capacity,
                                  bucket.size = bucket.size,
                                  use.smooth = use.smooth))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    print(paste0(i,' / ',iterations))
  }
  return(chain)
}

source('r/plot.mcmc.r')

# packages
library(doBy)
library(zoo)

# fit.mcmc.pace.func <- function(df = gcc.met.pace.df,
#                                species.in = 'Luc',prep.in = 'Control', 
#                                temp.in ='Ambient',subplot =NA,
#                                my.fun = phenoGrass.func.v11,
#                                out.nm.note = '',use.smooth = FALSE,
#                                day.lag = 3,
#                                bucket.size = 300,
#                                swc.capacity = 0.13,
#                                swc.wilt = 0.05,
#                                n.iter = 10000){
#   s.time <- Sys.time()
#   gcc.met.pace.df.16 <- get.pace.func(df,
#                                       species.in =species.in,
#                                       prep.in = prep.in,
#                                       temp.in =temp.in,
#                                       subplot = subplot)
#   gcc.met.pace.df.16$map <- 760
#   
#   # para values####
#   par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
#     f.t.opt = c(10,15,20,NA,NA,NA),
#     f.extract = c(0.05,0.075,0.1,NA,NA,NA),
#     f.sec = c(0.05,0.1,0.15,NA,NA,NA),
#     f.growth = c(0.1,0.2,0.3,NA,NA,NA),
#     q = c(0.5,1,2,NA,NA,NA))
#   row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
# 
#   # this assume 100% of the data falls into the max min range
#   # in a normal distribution for proposal.func
#   par.df['stdv',] <- ((par.df['max',] - par.df['min',])/100)
# 
#   # start mcmc fiting######
# 
#   # soil.water.var <- quantile(gcc.met.pace.df.16$vwc,c(.1,.99))
#   chain.fes=list()
#   for(n.chain in 1:3){
#     chain.fes[[n.chain]] = mh.MCMC.func(n.iter,
#                                         par.df,
#                                         gcc.met.pace.df.16,
#                                         bucket.size = bucket.size,
#                                         day.lay = day.lag,
#                                         swc.capacity = swc.capacity,
#                                         swc.wilt = swc.wilt,
#                                         my.fun = my.fun,
#                                         use.smooth = use.smooth)
#   }
# 
#   if(use.smooth==TRUE){
#     smooth.nm='sm'
#   }else{
#     smooth.nm=''
#   }
# 
#   if(is.na(subplot)){
#     out.name <- sprintf('cache/%s%schain.%s.%s.%s.rds',smooth.nm,out.nm.note,species.in,prep.in,temp.in)
#   }else{
#     out.name <- sprintf('cache/%schain.%s.rds',out.nm.note,subplot)
#   }
# 
#   saveRDS(chain.fes,out.name)
# 
#   print(Sys.time() - s.time)
# }

# wrpped function to fit mcmc to a df
fit.mcmc.2q.func <- function(df = gcc.met.pace.df,
                             species.in = 'Luc',prep.in = 'Control', 
                             temp.in ='Ambient',subplot =NA,
                             my.fun = phenoGrass.func.v13,
                             out.nm.note = '',use.smooth = FALSE,
                             day.lag = 3,
                             bucket.size = 300,
                             swc.capacity = 0.13,
                             swc.wilt = 0.05,
                             n.iter = 10000,
                             norm.min.max=NULL,
                             cal.initial=F,
                             par.df,q.given =NULL,q.s.given=NULL){
  
  # check for error in input
  if(ncol(par.df)==6 & !is.null(q.s.given)){
    stop('qs not set correctly. ')
  }

  if(ncol(par.df)==5 & !is.null(q.given)){
    stop('q not set correctly. ')
  }
  
  # 
  s.time <- Sys.time()
  gcc.met.pace.df.16 <<- get.pace.func(df,
                                      species.in =species.in,
                                      prep.in = prep.in,
                                      temp.in =temp.in,
                                      subplot = subplot,
                                      norm.min.max = norm.min.max)
  
  my.fun <<-  my.fun
  use.smooth <<-  use.smooth
  day.lag <<- day.lag
  bucket.size <<-  bucket.size
  swc.capacity <<- swc.capacity
  swc.wilt <<-  swc.wilt
  n.iter <<-  n.iter
  
  if(cal.initial){
    source('r/deoptimal_initial.R')
    initial.vec <- get.ini.func(par.df = par.df,q.given =q.given,q.s.given=q.s.given)
     out.nm <- paste0('tmp/deopt_',species.in,prep.in,temp.in,'.rds')
     saveRDS(initial.vec,out.nm)
    }else{
    initial.vec<-NULL
  }

  # this assume 100% of the data falls into the max min range
  # in a normal distribution for proposal.func

  if(is.null(initial.vec)){
    par.df['initial',] <- c(20,0.5,0.005,0.15,3,0.5)[seq_along(par.df['initial',])]
  }else{
    par.df['initial',] <- initial.vec
    # par.df['max',] <- initial.vec*2
    # par.df['min',] <- initial.vec/2
  }
  par.df['stdv',] <- (par.df['max',] - par.df['min',])/100
  
  par.df <<- par.df
  # print(par.df)
  
  # gcc.met.pace.df.16$map <- 760
  
  # # para values####
  # par.df <<- data.frame(#f.h = c(200,220,240,NA,NA),
  #   f.t.opt = c(10,15,20,NA,NA,NA),
  #   f.extract = c(0.05,0.075,0.1,NA,NA,NA),
  #   f.sec = c(0.05,0.1,0.15,NA,NA,NA),
  #   f.growth = c(0.1,0.2,0.3,NA,NA,NA),
  #   q = c(0.5,1,2,NA,NA,NA),
  #   q.s = c(0.5,1,2,NA,NA,NA))
  # row.names(par.df) <<- c('min','initial','max','fit','stdv','prop')
  # 
  # # this assume 100% of the data falls into the max min range
  # # in a normal distribution for proposal.func
  # par.df['stdv',] <<- (par.df['max',] - par.df['min',])/100
  
  # # start mcmc fiting######
  # # soil.water.var <- quantile(gcc.met.pace.df.16$vwc,c(.1,.99))
  # chain.fes=list()
  # for(n.chain in 1:3){
  #   chain.fes[[n.chain]] = mh.MCMC.func.2q(n.iter,
  #                                          par.df,
  #                                          gcc.met.pace.df.16,
  #                                          bucket.size = bucket.size,
  #                                          day.lay = day.lag,
  #                                          swc.capacity = swc.capacity,
  #                                          swc.wilt = swc.wilt,
  #                                          my.fun = my.fun,
  #                                          use.smooth = use.smooth)
  # }

  
  #setup parallel backend to use many processors
  cores=3#detectCores(logical = FALSE)
  cl <- makeCluster(cores[1]) #not to overload your computer
  registerDoParallel(cl)
  #stop cluster
  on.exit(stopCluster(cl))
  
  chain.fes <- foreach(i=1:3,
                       .packages = 'mvtnorm',
                       .export=ls(envir=globalenv())) %dopar% {
                         
                         chain.tmp = mh.MCMC.func.2q(iterations=n.iter,
                                                     par.df = par.df,
                                                     gcc.met.pace.df.16 = gcc.met.pace.df.16,
                                                     bucket.size = bucket.size,
                                                     day.lay = day.lag,
                                                     swc.capacity = swc.capacity,
                                                     swc.wilt = swc.wilt,
                                                     my.fun = phenoGrass.func.v13,
                                                     use.smooth = T,q.given =q.given,q.s.given=q.s.given)
                       }

  # stopCluster(cl)
  
  # save file
  if(use.smooth==TRUE){
    smooth.nm='sm'
  }else{
    smooth.nm=''
  }
  
  if(is.na(subplot)){
    out.name <- sprintf('cache/%s%schain.%s.%s.%s.rds',smooth.nm,out.nm.note,species.in,prep.in,temp.in)
  }else{
    out.name <- sprintf('cache/%schain.%s.rds',out.nm.note,subplot)
  }
  
  saveRDS(chain.fes,out.name)
  
  time.used <- (Sys.time() - s.time)
  
  df <- data.frame(site = species.in,
                   time = time.used,
                   when = Sys.time())
  
  write.table(df, file = "time.used.txt", sep = "\t",append = T,
              row.names = FALSE)
}

# generic function to do mcmc in paralle
mh.MCMC.func.2q <- function(iterations,par.df,
                            gcc.met.pace.df.16,
                            bucket.size,swc.wilt ,
                            swc.capacity,day.lay,
                            my.fun,
                            use.smooth,q.given,q.s.given){


  
  # get prior 
  prior.prob <- prior.func(par.df)
  
  # start MC ####
  # intial
  chain = array(dim = c(iterations+1,ncol(par.df)))
  # chain[1,] = as.numeric(par.df['initial',])
  chain[1,] = proposal.func(par.df['initial',],par.df)
  
  # chain move on
  for (i in 1:iterations){
    # prpose a set of par values based on previous chain value
    proposal = proposal.func(chain[i,],par.df)
    
    # deal with constant sensitivity
    # if(length(proposal)<6){
    #   q.s.val = 0
    # }else{
    #   q.s.val = proposal[6]
    # }
    # 
    # if(length(proposal)<5){
    #   q.val = 0
    # }else{
    #   q.val = proposal[5]
    # }
    
    if(is.null(q.given)){
      q.val = proposal[5]
    }else{
      q.val = q.given
    }
    
    if(is.null(q.s.given)){
      q.s.val = proposal[6]
    }else{
      q.s.val = q.given
    }

    # prior.prob,data,data.sd,bucket.size = 300,...
    probab = exp(posterior.func(prior.prob,FUN = my.fun,
                                gcc.df = gcc.met.pace.df.16,
                                f.h = 222,
                                f.t.opt = proposal[1],
                                f.extract = proposal[2],
                                f.sec = proposal[3],
                                f.growth = proposal[4] ,
                                q = q.val ,
                                q.s = q.s.val ,
                                t.max = 45,
                                day.lay = day.lay,
                                bucket.size = bucket.size,
                                swc.wilt = swc.wilt ,
                                swc.capacity = swc.capacity,
                                use.smooth = use.smooth) - 
                   posterior.func(prior.prob,FUN = my.fun,
                                  gcc.df = gcc.met.pace.df.16,
                                  f.h = 222,
                                  f.t.opt = chain[i,1],
                                  f.extract = chain[i,2],
                                  f.sec = chain[i,3],
                                  f.growth = chain[i,4] ,
                                  q = q.val ,
                                  q.s = q.s.val ,
                                  t.max = 45,
                                  day.lay = day.lay,
                                  swc.wilt = swc.wilt ,
                                  swc.capacity = swc.capacity,
                                  bucket.size = bucket.size,
                                  use.smooth = use.smooth))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    # print(paste0(i,' / ',iterations))
  }
  return(chain)
}
