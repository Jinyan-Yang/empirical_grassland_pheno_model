source('r/read_spc_nm.R')
# process ym data and pace data
param.ls <- list()

for (i in seq_along(species.vec)) {
  fn.con.11 <- sprintf('cache/smv1.2q.chain.%s.Control.Ambient.rds',
                       species.vec[i])
  
  par.con.11 <- readRDS(fn.con.11)
  
  chain.3.ls.new = lapply(par.con.11,function(m.in)m.in[round((0.75)*nrow(m.in)):nrow(m.in),])
  
  chain.fes <- do.call(rbind,chain.3.ls.new)
  
  # chain.fes$spc <- species.vec.nm[i]
  
  param.ls[[i]] <-  subset(chain.fes,select = -c(ll))
  
}

names(param.ls) <- species.vec.nm

out.ls <- lapply(param.ls,cor)

# cov###
for (i in seq_along(out.ls)) {
  out.ls[[i]] <- rbind(rep(NA,6),out.ls[[i]]
                       )
  
  rownames(out.ls[[i]])[1] <- names(out.ls)[i]
}

out.df <- do.call(rbind,out.ls)

write.csv(out.df,'cache/cor_matrix.csv')
# 
png('figures/cov_plot.png',width = 800,height = 800*.618)
par(mfrow=c(2,2))
# plot###
library(vioplot)
palette(c(col.df$iris))
col.nm.vec <- c(1,1,2,2,2,3,3,3,4,4)
pch.vec <- c(16,3,16,3,2,16,3,2,16,3)

# senesce
plot(f.sec~q.s,data = param.ls[[1]],pch=pch.vec[1],col=col.nm.vec[1],xlim=c(0,15),ylim=c(0.01,0.5),
       ylab = expression(r[senescence]),xlab=expression(q[senescence]))
 for (i in 2:10) { 
  points(f.sec~q.s,data = param.ls[[i]],pch=pch.vec[i],col=col.nm.vec[i])
}

legend('bottomright',legend = species.vec.nm,pch = pch.vec,col=col.nm.vec,ncol=2)
legend('topleft',legend = '(a)',bty='n')

# growth
plot(f.growth~q,data = param.ls[[1]],pch=pch.vec[1],col=col.nm.vec[1],xlim=c(0,15),ylim=c(0.01,0.5),
     ylab = expression(r[growth]),xlab=expression(q[growth]))
for (i in 2:10) { 
  points(f.growth~q,data = param.ls[[i]],pch=pch.vec[i],col=col.nm.vec[i])
}

# legend('bottomright',legend = species.vec.nm,pch = pch.vec,col=col.nm.vec,ncol=2)
legend('topleft',legend = '(b)',bty='n')

# q growth vs sene
plot(q.s~q,data = param.ls[[1]],pch=pch.vec[1],col=col.nm.vec[1],xlim=c(0,15),ylim=c(0,15),
     ylab = expression(q[senescence]),xlab=expression(q[growth]))
for (i in 2:10) { 
  points(q.s~q,data = param.ls[[i]],pch=pch.vec[i],col=col.nm.vec[i])
}

# legend('bottomright',legend = species.vec.nm,pch = pch.vec,col=col.nm.vec,ncol=2)
legend('topleft',legend = '(c)',bty='n')

# growth vs senc
plot(f.growth~f.sec,data = param.ls[[1]],pch=pch.vec[1],col=col.nm.vec[1],xlim=c(0,0.5),ylim=c(0,0.5),
     ylab = expression(r[growth]),xlab=expression(r[senescence]))
for (i in 2:10) { 
  points(f.growth~f.sec,data = param.ls[[i]],pch=pch.vec[i],col=col.nm.vec[i])
}

# legend('bottomright',legend = species.vec.nm,pch = pch.vec,col=col.nm.vec,ncol=2)
legend('topleft',legend = '(d)',bty='n')

dev.off()
  