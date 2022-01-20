# 
source('r/load.R')
source('r/read_spc_nm.R')

# read fitted chain the de debug folder
fn <- 'debug/smtest.flux.chain.flux.Control.Ambient.rds'

in.chain =  readRDS(fn)
# 
chain.fes <- do.call(rbind,in.chain)
chain.fes <- chain.fes[!duplicated(chain.fes),]

#read met 
ym.18.df <- get.ym.func(18)

gcc.met.con.df <- get.paddock.func('control')
#need to specif i for code to work####
# ym i=9; flux tower, i=10

# raad met data for site
if(species.vec[i]=='ym'){
  df = ym.18.df
  # 
  # c.wd <- getwd()
  # setwd('c:/repo/dn_gcc/')
  ym.met.df <- readRDS('cache/ym/ym.met.rds')
  # 
  swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
  swc.cap = round(swc.ym.con[[2]]*10)/10
  swc.wilt = round(swc.ym.con[[1]]*100)/100
  bucket.size=1000
}else if(species.vec[i]=='flux'){
  df = gcc.met.con.df
  swc.q.con <- quantile(gcc.met.con.df$vwc,na.rm=T,probs = c(0.01,0.99))
  swc.cap =  round(swc.q.con[[2]]*10)/10
  swc.wilt = round(swc.q.con[[1]]*100)/100
  bucket.size=1000
}else{
  df = gcc.met.pace.df
  swc.cap = 0.13
  swc.wilt = 0.05
  bucket.size=300
}

# make plots
pdf('test10k.pdf',width = 8,height = 8*.618)
# plot time series
plot.mcmc.func.2q(df,species.vec[i],
                  prep.in='Control',temp.in='Ambient',
                  my.fun = phenoGrass.func.v13,
                  nm.note='smv13.2q.',use.smooth = TRUE,day.lag = 3,
                  swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,
                  bucket.size = bucket.size)

# plot chains with given num of iteration
for(par.i in 1:7){
  plot.line.mcmc.func(in.chain,val.nm=par.i,range.iter = 17000:20000)
}

# histgram
plot.check.mcmc.func(chain.fes)


dev.off()
