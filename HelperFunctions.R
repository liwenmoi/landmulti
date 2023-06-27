# some auxiliary functions
s1 <- function(x) x
s2 <- function(x) x
delta <- function(x) x

# logit link function
g.logit <- function(xx){exp(xx)/(1+exp(xx))}

# calculate kernel matrix
Kern.FUN.log  <- function(S1i,ti,hi) {
  out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
  norm.k = dnorm(out)/hi
  norm.k
}

# calculate kernel matrix using bivariate normal
Kern.FUN.2.log  <- function(S1i,S2i, ti,Hi) ## returns an (n x nz) matrix
{ out = (log(rbind(S1i,S2i))- log(ti))
norm.k =emdbook::dmvnorm(x=t(out), mu=c(0,0), Sigma=Hi)
norm.k = det(Hi)^(-1/2) * norm.k
norm.k
}

################Function 3 bandwidth for group2#################################
min.BW.cens.ex.gr2 <- function(da2,t0,L,h,folds,reps,s.seq)
{
Y     = da2$Y
XS1   = da2$XS1
XS2   = da2$XS2
delta = da2$delta
X     = data.matrix(subset(da2, select = -c(Y, XS1, XS2, delta, what, group)))
Wi    = da2$what

# logit link function
g.logit <- function(xx){exp(xx)/(1+exp(xx))}

# calculate kernel matrix (a log is added. Need to discuss)
Kern.FUN.log  <- function(S1i,ti,hi) {
  out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
  norm.k = dnorm(out)/hi
  norm.k
}
loc.fun.ex <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS1 <=t0)
  K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}
n = dim(da2)[1]; nv = floor(n/folds)
BW.vec = vector(length=folds)
replic = matrix(nrow = reps, ncol =folds)
for(lay in 1:reps) {
  for(k in 1:folds) {
    ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
    # subset.val = ind.val[dS1[ind.val]==1 & XS1[ind.val]<= t0 & Y[ind.val] > t0];
    # subset.tra = ind.tra[dS1[ind.tra]==1 & XS1[ind.tra]<= t0 & Y[ind.tra] > t0];
    # calculate the est.coef for P in the MSE
    P.md= loc.fun.ex(s.seq = s.seq[1], data=da2[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
    p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
    BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(Wi[ind.val]))
    BW.vec[k] = BW
  }
  replic[lay,] = BW.vec
}
mean(replic)
}
loc.fun.ex <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS1 <=t0)
  K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}
################Function 3 bandwidth for group3#################################
min.BW.cens.ex.gr3 <- function(da3,t0,L,h,folds,reps, s.seq)
{
  Y     = da3$Y
  XS1   = da3$XS1
  XS2   = da3$XS2
  delta = da3$delta
  X     = data.matrix(subset(da3, select = -c(Y, XS1, XS2, delta, what, group)))
  Wi    = da3$what


  # logit link function
  g.logit <- function(xx){exp(xx)/(1+exp(xx))}

  # calculate kernel matrix (a log is added. Need to discuss)
  Kern.FUN.log  <- function(S1i,ti,hi) {
    out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
    norm.k = dnorm(out)/hi
    norm.k
  }

  loc.fun.ex.gr3 <- function(s.seq, data, t0,  L, h, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
    Wi    = data$what

    fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

    # `index.sub` can be removed, but it does not do any harm to leave it here...
    index.sub = (data$Y > t0) & (data$XS2 <= t0)
    K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
    est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(s.seq)) {
      m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da3)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
      # subset.val = ind.val[dS2[ind.val]==1 & XS2[ind.val]<= t0 & Y[ind.val] > t0];
      # subset.tra = ind.tra[dS2[ind.tra]==1 & XS2[ind.tra]<= t0 & Y[ind.tra] > t0];
      # s.seq = s.seq
      # calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr3(s.seq = s.seq[1], data=da3[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}
loc.fun.ex.gr3 <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS2 <= t0)
  K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}
################Function 3 bandwidth for group4#################################
min.BW.cens.ex.gr4 <- function(da4,t0,L,h,folds,reps, s.seq)
{
  H=rbind(c(h[1],0), c(0,h[2]))
  Y     = da4$Y
  XS1   = da4$XS1
  XS2   = da4$XS2
  delta = da4$delta
  X     = data.matrix(subset(da4, select = -c(Y, XS1, XS2, delta, what, group)))
  Wi    = da4$what

  g.logit <- function(xx){exp(xx)/(1+exp(xx))}
  Kern.FUN.2.log  <- function(S1i,S2i, ti,Hi) ## returns an (n x nz) matrix
  { out = (log(rbind(S1i,S2i))- log(ti))
  norm.k =emdbook::dmvnorm(x=t(out), mu=c(0,0), Sigma=Hi)
  norm.k = det(Hi)^(-1/2) * norm.k
  norm.k
  }

  loc.fun.ex.gr4 <- function(s.seq, data, t0,  L, H, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    dS1   = data$deltaS1
    dS2   = data$deltaS2
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
    Wi    = data$what

    fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
    #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
    K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
    K=t(as.matrix(K))
    est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(t(as.matrix(s.seq))[,1])) {
      m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }


  n = dim(da4)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)

      #calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr4(s.seq = s.seq, data=da4[ind.tra,], t0=t0, L=L, H=H, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}
loc.fun.ex.gr4 <- function(s.seq, data, t0,  L, H, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  dS1   = data$deltaS1
  dS2   = data$deltaS2
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group))
  Wi    = data$what

  fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
  #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
  K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
  K=t(as.matrix(K))
  est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(t(as.matrix(s.seq))[,1])) {
    m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}




###############
min.BW.cens.ex.gr2.SE <- function(da2,t0,L,h,folds,reps,s.seq)
{
  Y     = da2$Y
  XS1   = da2$XS1
  XS2   = da2$XS2
  delta = da2$delta
  X     = data.matrix(subset(da2, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap)))
  Wi    = da2$what

  # logit link function
  g.logit <- function(xx){exp(xx)/(1+exp(xx))}

  # calculate kernel matrix (a log is added. Need to discuss)
  Kern.FUN.log  <- function(S1i,ti,hi) {
    out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
    norm.k = dnorm(out)/hi
    norm.k
  }

  loc.fun.ex.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
    Wi    = data$what

    fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

    # `index.sub` can be removed, but it does not do any harm to leave it here...
    index.sub = (data$Y > t0) & (data$XS1 <=t0)
    K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
    est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(s.seq)) {
      m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da2)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
      # subset.val = ind.val[dS1[ind.val]==1 & XS1[ind.val]<= t0 & Y[ind.val] > t0];
      # subset.tra = ind.tra[dS1[ind.tra]==1 & XS1[ind.tra]<= t0 & Y[ind.tra] > t0];
      # calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.SE(s.seq = s.seq[1], data=da2[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}
min.BW.cens.ex.gr3.SE <- function(da3,t0,L,h,folds,reps, s.seq)
{
  Y     = da3$Y
  XS1   = da3$XS1
  XS2   = da3$XS2
  delta = da3$delta
  X     = data.matrix(subset(da3, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap)))
  Wi    = da3$what

  # logit link function
  g.logit <- function(xx){exp(xx)/(1+exp(xx))}

  # calculate kernel matrix (a log is added. Need to discuss)

  Kern.FUN.log  <- function(S1i,ti,hi) {
    out = (log(landpred::VTM(vc =S1i, dm=length(ti)))- log(ti))/hi
    norm.k = dnorm(out)/hi
    norm.k
  }

  loc.fun.ex.gr3.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
    Wi    = data$what

    fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

    # `index.sub` can be removed, but it does not do any harm to leave it here...
    index.sub = (data$Y > t0) & (data$XS2 <= t0)
    K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
    est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(s.seq)) {
      m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da3)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)
      # subset.val = ind.val[dS2[ind.val]==1 & XS2[ind.val]<= t0 & Y[ind.val] > t0];
      # subset.tra = ind.tra[dS2[ind.tra]==1 & XS2[ind.tra]<= t0 & Y[ind.tra] > t0];
      # s.seq = s.seq
      # calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr3.SE(s.seq = s.seq[1], data=da3[ind.tra,], t0=t0, L=L, h=h, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}
min.BW.cens.ex.gr4.SE <- function(da4,t0,L,h,folds,reps, s.seq)
{
  H=rbind(c(h[1],0), c(0,h[2]))
  Y     = da4$Y
  XS1   = da4$XS1
  XS2   = da4$XS2
  delta = da4$delta
  X     = data.matrix(subset(da4, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap)))
  Wi    = da4$what

  g.logit <- function(xx){exp(xx)/(1+exp(xx))}
  Kern.FUN.2.log  <- function(S1i,S2i, ti,Hi) ## returns an (n x nz) matrix
  { out = (log(rbind(S1i,S2i))- log(ti))
  norm.k =emdbook::dmvnorm(x=t(out), mu=c(0,0), Sigma=Hi)
  norm.k = det(Hi)^(-1/2) * norm.k
  norm.k
  }
  loc.fun.ex.gr4.SE <- function(s.seq, data, t0,  L, H, weight = NULL) {
    Y     = data$Y
    XS1   = data$XS1
    XS2   = data$XS2
    delta = data$delta
    dS1   = data$deltaS1
    dS2   = data$deltaS2
    X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
    Wi    = data$what

    fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
    #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
    K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
    K=t(as.matrix(K))
    est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
    invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
    con.mat=matrix(nrow=length(s.seq), ncol = 1)
    for(i in 1:length(t(as.matrix(s.seq))[,1])) {
      m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
      con.mat[i,]=m$converged
    }
    return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
  }

  n = dim(da4)[1]; nv = floor(n/folds)
  BW.vec = vector(length=folds)
  replic = matrix(nrow = reps, ncol =folds)
  for(lay in 1:reps) {
    for(k in 1:folds) {
      ind.val = sample(1:n, nv); ind.tra = setdiff(1:n,ind.val)

      #calculate the est.coef for P in the MSE
      P.md= loc.fun.ex.gr4.SE(s.seq = s.seq, data=da4[ind.tra,], t0=t0, L=L, H=H, weight = Wi[ind.tra])
      p.hat = g.logit(P.md$est.mat[,1] + P.md$est.mat[,-1]%*%t(X[ind.val,]))
      BW = sum((((Y[ind.val]<=t0+L) - p.hat)^2)*(Y[ind.val]>t0)*(XS1[ind.val]<=(t0))*(XS2[ind.val]<=(t0))* (Wi[ind.val]))
      BW.vec[k] = BW
    }
    replic[lay,] = BW.vec
  }
  mean(replic)
}
loc.fun.ex.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
  Wi    = data$what

  fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS1 <=t0)
  K = Kern.FUN.log(S1i=XS1,ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml2, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}
loc.fun.ex.gr3.SE <- function(s.seq, data, t0,  L, h, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
  Wi    = data$what

  fml3 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(names(X), collapse = "+")))

  # `index.sub` can be removed, but it does not do any harm to leave it here...
  index.sub = (data$Y > t0) & (data$XS2 <= t0)
  K = Kern.FUN.log(S1i=XS2[index.sub],ti=s.seq,hi=h)
  est.mat = matrix(nrow=length(s.seq), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(s.seq), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(s.seq)) {
    m = glm(fml3, data = data, family = "binomial", weights = K[i,]*weight[index.sub])
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}
loc.fun.ex.gr4.SE <- function(s.seq, data, t0,  L, H, weight = NULL) {
  Y     = data$Y
  XS1   = data$XS1
  XS2   = data$XS2
  delta = data$delta
  dS1   = data$deltaS1
  dS2   = data$deltaS2
  X     = subset(data, select = -c(Y, XS1, XS2, delta, what, group, Expwt, W.resap))
  Wi    = data$what

  fml4 <- as.formula(paste("1*(data$Y<= t0+L) ~ ", paste(names(X), collapse = "+")))
  #index.sub = (data$Y > t0) & (data$XS1 <= t0) & (data$deltaS1 == 1) & (data$XS2 <= t0) & (data$deltaS2 == 1)
  K = Kern.FUN.2.log(S1i=XS1, S2i=XS2, ti=s.seq, Hi=H)
  K=t(as.matrix(K))
  est.mat    = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = dim(as.matrix(X))[2] + 1)
  invinf.mat = matrix(nrow=length(t(as.matrix(s.seq))[,1]), ncol = (dim(as.matrix(X))[2] + 1)^2)
  con.mat=matrix(nrow=length(s.seq), ncol = 1)
  for(i in 1:length(t(as.matrix(s.seq))[,1])) {
    m = glm(fml4, data = data, family = "binomial", weights = K[i,]*weight)
    est.mat[i,] = m$coeff
    invinf.mat[i,] = as.vector(vcov(m))
    con.mat[i,]=m$converged
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat, "conv"=con.mat))
}





multipred <- function(data, formula, t0, L,
                      SE=FALSE, SE.gs = FALSE,
                      s1_beta1=NULL, s2_beta2=NULL, s1s2_beta3=NULL,
                      grid1 = seq(0.01, 5, length.out=20),
                      grid2 = seq(0.01, 5, length.out=20),
                      grid3 = list(seq(0.01, 5, length.out=20),
                                   seq(0.01, 5, length.out=20)),
                      folds.grid = 8, reps.grid = 3,
                      c01 = 0.1,
                      c02 = 0.1,
                      c03 = 0.05,
                      B = 500,
                      gs.method = "loop",
                      gs.cl = NULL
) {

  # ------------------ setups for parallel computing -------------------------#
  system <- Sys.info()['sysname']
  cl <- NULL
  if (system == 'Windows') {
    cl <- makeCluster(getOption('cl.cores', gs.cl))
    registerDoParallel(cl)
    registerDoSEQ()
    on.exit(stopCluster(cl))
  } else {
    options('mc.cores' = cores)
    registerDoParallel(cores)
  }

  s1 <- function(x) x
  s2 <- function(x) x
  delta <- function(x) x

  # -------------------------------------------------------#
  # Qian's new code for data manipulation
  model.frame(formula, data)
  mf <- model.frame(formula, data)

  pos_s1 <- grep("s1", names(mf))
  pos_s2 <- grep("s2", names(mf))
  # pos_delta <- grep("delta", names(mf))

  XS1 <- mf[[pos_s1]]
  XS2 <- mf[[pos_s2]]

  delta <- mf[[1]][,2]

  Y <- mf[[1]][,1]
  X <- mf[-c(1, pos_s1, pos_s2)]

  Xnames <- names(mf)[-c(1, pos_s1, pos_s2)]

  mydata <- as.data.frame(cbind(Y, delta, XS1, XS2, X))

  mydata$what=Wi.FUN(data=mydata, t0=t0, tau=L, weight.given=NULL)
  mydata$group=0
  mydata$group[(mydata$Y>t0)&(mydata$XS1>t0)&(mydata$XS2>t0)]=1
  mydata$group[(mydata$Y>t0)&(mydata$XS1<=t0)&(mydata$XS2>t0)]=2
  mydata$group[(mydata$Y>t0)&(mydata$XS1>t0)&(mydata$XS2<=t0)]=3
  mydata$group[(mydata$Y>t0)&(mydata$XS1<=t0)&(mydata$XS2<=t0)]=4

  da1 <- subset(mydata, group==1)
  da2 <- subset(mydata, group==2)
  da3 <- subset(mydata, group==3)
  da4 <- subset(mydata, group==4)
  n.f.b=nrow(mydata)

  #-------------------- calculate beta_t0 for group 1 -------------------------#

  fml1 <- as.formula(paste("1*(Y <= t0+L) ~", paste(Xnames, collapse = "+")))
  G1_model = glm(fml1,  data=da1, family = "binomial", weights = what)
  G1_coef = G1_model$coeff
  G1_con= G1_model$conv
  res1 <- list("group1_coef" = G1_coef,
               "group1_convergence" = G1_con)

  res <- list("group1" = res1)

  if (!is.null(s1_beta1)) {

    #-------------------- calculate beta_s1 for group2 ---------------------------#
    # find the bandwidth that optimize the MSE
    fml2 <- as.formula(paste("1*(data$Y[index.sub]<= t0+L)~ ", paste(Xnames, collapse = "+")))

    h2.min1 = gridSearch(fun = min.BW.cens.ex.gr2, da2=da2,t0=t0,L=L, s.seq = s1_beta1,
                         folds= folds.grid, reps= reps.grid, npar = 1, levels = list(a=grid1),
                         method = gs.method, cl=gs.cl)

    # shrink bandwidth to get the final bandwidth
    h2.final1= h2.min1$minlevels/(n.f.b^c01)

    # apply final bandwidth to the estimating equation
    G2_model1 = loc.fun.ex(s.seq=s1_beta1, data=da2, t0=t0,  L=L, h=h2.final1, weight = da2$what)
    G2_coef1 = G2_model1$est.mat
    G2_con1= G2_model1$conv


    res2 <- list("group2_coef" = G2_coef1,
                 "group2_convergence" = G2_con1)
    colnames(res2$group2_coef) <- c("Intercept", Xnames)
    res <- append(res, res2)
  }

  #-------------------- calculate beta_s2 for group3 ---------------------------#
  # find the bandwidth that optimize the MSE

  if (!is.null(s2_beta2)) {

    h3.min1 = gridSearch(fun = min.BW.cens.ex.gr3, da3=da3,t0=t0, L=L,
                         folds= folds.grid, reps= reps.grid, npar = 1, s.seq=s2_beta2, levels = list(a=grid2),
                         method = gs.method, cl=gs.cl)

    # shrink bandwidth to get the final bandwidth
    h3.final1= h3.min1$minlevels/(n.f.b^c02)

    # apply final bandwidth to the estimating equation
    G3_model1 = loc.fun.ex.gr3(s.seq=s2_beta2, data=da3, t0=t0,  L=L, h=h3.final1, weight = da3$what)
    G3_coef1 = G3_model1$est.mat
    G3_con1= G3_model1$conv

    res3 <- list("group3_coef" = G3_coef1,
                 "group3_convergence" = G3_con1)
    colnames(res3$group3_coef) <- c("Intercept", Xnames)
    res <- append(res, res3)
  }

  #-------------------- calculate beta_s3 for group4 ---------------------------#
  # find the bandwidth that optimize the MSE

  if (!is.null(s1s2_beta3)) {
    res4 <- c()
    res4c <- c()

    # convert s1s2_beta3 into required format
    if (class(s1s2_beta3)[1] == "numeric") {
      s1s2_beta3 <- matrix(s1s2_beta3, nrow = 1)} else {
        s1s2_beta3 <- as.matrix(s1s2_beta3)
      }

    for (i in 1:nrow(s1s2_beta3)) {
      status4_11=gridSearch(fun=min.BW.cens.ex.gr4,da4=da4,t0=t0,L=L,
                            s.seq = s1s2_beta3[i,],
                            npar= 2, folds=folds.grid, reps=reps.grid,
                            levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
                            method = gs.method, cl=gs.cl)
      h4.band11=status4_11$minlevels
      h4.final11=(sqrt(h4.band11)/(n.f.b^c03))^2
      G4_model11 = loc.fun.ex.gr4(s.seq=s1s2_beta3[i,], data=da4, t0=t0,
                                  L=L, H=rbind(c(h4.final11[1],0), c(0,h4.final11[2])), weight = da4$what)
      G4_coef11 = G4_model11$est.mat
      G4_con11 = G4_model11$conv[1,1]
      res4i <- G4_coef11
      res4ic <- G4_con11
      res4 <- rbind(res4, res4i)
      res4c <-rbind(res4c, res4ic)
      # colnames(res4$group4_coef) <- c("Intercept", Xnames)
    }
    res4 <- list("group4_coef" = res4, "group4_convergence" = res4c)
    colnames(res4$group4_coef) <- c("Intercept", Xnames)
    res <- append(res, res4)

  }

  #------------ Calculate empirical variance if `SE == T` ---------------------#
  if(SE) {

    if(SE.gs == T){
      resap.SE.vec <- c()
      for (i in 1:B) {
        mydata.o=mydata
        mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
        mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
        mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
        mydata.o$group=0
        mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2>t0)]=1
        da1.sap=mydata.o[mydata.o$group==1,]

        G1_model.sap = glm(fml1,  data=da1.sap, family = "binomial", weights = da1.sap$W.resap)
        G1_coef.sap = G1_model.sap$coeff
        resap.SE.vec=rbind(resap.SE.vec, G1_coef.sap)
      }
      resap.SE.vec=as.data.frame(resap.SE.vec)
      group1_SE=apply(resap.SE.vec, 2, sd)
      group1_SE <- as.data.frame(group1_SE)
      res <- append(res, group1_SE)
      names(res$group1_SE) = c("Intercept", Xnames)

      #------------------ SE for group2 ------------------#
      if (!is.null(s1_beta1)) {
        se.list2 <- list()
        for (i in 1:length(s1_beta1)) {
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2>t0)]=2
            da2.sap=mydata.o[mydata.o$group==2,]

            h2.min1 = gridSearch(fun = min.BW.cens.ex.gr2.SE, da2=da2.sap,t0=t0, L=L,
                                 folds= folds.grid, reps= reps.grid, npar = 1,
                                 s.seq=s1_beta1, levels = list(a=grid1),
                                 method = gs.method, cl=gs.cl)
            # shrink bandwidth to get the final bandwidth
            h2.final1.sap= h2.min1$minlevels/(nrow(da2.sap)^c02)
            G2_model1.sap = loc.fun.ex.SE(s.seq=s1_beta1[i], data=da2.sap, t0=t0,  L=L, h=h2.final1.sap,
                                          weight = da2.sap$W.resap)
            G2_coef1.sap = G2_model1.sap$est.mat
            resap.SE.vec=rbind(resap.SE.vec, G2_coef1.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group2.row",i)
          se.list2[[df_name]] <- SE

          #resap.SE.seq = append(resap.SE.seq, G2_coef1.sap)
        }
        res <- append(res, se.list2)
      }
      #------------------ SE for group3 ------------------#
      if (!is.null(s2_beta2)) {

        se.list3 <- list()

        for (i in 1:length(s2_beta2)){
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2<=t0)]=3
            da3.sap=mydata.o[mydata.o$group==3,]

            h3.min1 = gridSearch(fun = min.BW.cens.ex.gr3.SE, da3=da3.sap,t0=t0, L=L,
                                 folds= folds.grid, reps= reps.grid, npar = 1,
                                 s.seq=s2_beta2, levels = list(a=grid2),
                                 method = gs.method, cl=gs.cl)
            # shrink bandwidth to get the final bandwidth
            h3.final1.sap= h3.min1$minlevels/(nrow(da3.sap)^c02)


            G3_model1.sap=loc.fun.ex.gr3.SE(s.seq=s2_beta2[i], data=da3.sap, t0=t0,  L=L, h=h3.final1.sap,
                                              weight = da3.sap$W.resap)
            G3_coef1.sap = G3_model1.sap$est.mat
            resap.SE.vec = rbind(resap.SE.vec, G3_coef1.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group3.row",i)
          se.list3[[df_name]] <- SE
        }
        res <- append(res, se.list3)
      }
      #------------------ SE for group4 ------------------#
      if (!is.null(s1s2_beta3)) {

        se.list4 <- list()

        for (i in 1:nrow(s1s2_beta3)) {
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2<=t0)]=4
            da4.sap=mydata.o[mydata.o$group==4,]

            status4_11=gridSearch(fun=min.BW.cens.ex.gr4.SE,da4=da4.sap,t0=t0,L=L,
                                  s.seq = s1s2_beta3[i,],
                                  npar= 2, folds=folds.grid, reps=reps.grid,
                                  levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
                                  method = gs.method, cl=gs.cl)
            h4.band11_sap=status4_11$minlevels
            h4.final11_sap=(sqrt(h4.band11_sap)/(nrow(da4.sap)^.05))^2
            G4_model11.sap = loc.fun.ex.gr4.SE(s.seq=as.matrix(s1s2_beta3)[1,], data=da4.sap, t0=t0,
                                               L=L, H=rbind(c(h4.final11_sap[1],0), c(0,h4.final11_sap[2])), weight = da4.sap$W.resap)
            G4_coef11.sap = G4_model11.sap$est.mat
            resap.SE.vec = rbind(resap.SE.vec, G4_coef11.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group4.row",i)
          se.list4[[df_name]] <- SE
        }
        res <- append(res, se.list4)
      }

    } else if(SE.gs == F)
      {
      se.list1 <- list()
      resap.SE.vec <- c()
      for (i in 1:B) {
        # print(i)
        mydata.o=mydata
        mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
        mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
        mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
        mydata.o$group=0
        mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2>t0)]=1
        da1.sap=mydata.o[mydata.o$group==1,]

        G1_model.sap = glm(fml1,  data=da1.sap, family = "binomial", weights = da1.sap$W.resap)
        G1_coef.sap = G1_model.sap$coeff
        resap.SE.vec=rbind(resap.SE.vec, G1_coef.sap)
      }
      resap.SE.vec=as.data.frame(resap.SE.vec)
      group1_SE=apply(resap.SE.vec, 2, sd)
      group1_SE <- as.data.frame(group1_SE)
      res <- append(res, group1_SE)
      names(res$group1_SE) = c("Intercept", Xnames)

    #------------------ SE for group2 ------------------#
      if (!is.null(s1_beta1)) {

        se.list2 <- list()

        for (i in 1:length(s1_beta1)) {
          resap.SE.vec <- c()
          for (j in 1:B) {
            mydata.o=mydata
            mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
            mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
            mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
            mydata.o$group=0
            mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2>t0)]=2
            da2.sap=mydata.o[mydata.o$group==2,]

            G2_model1.sap = loc.fun.ex.SE(s.seq=s1_beta1[i], data=da2.sap, t0=t0,  L=L, h=h2.final1,
                                        weight = da2.sap$W.resap)
            G2_coef1.sap = G2_model1.sap$est.mat
            resap.SE.vec=rbind(resap.SE.vec, G2_coef1.sap)
          }
          resap.SE.vec=as.data.frame(resap.SE.vec)
          SE=apply(resap.SE.vec, 2, sd)
          names(SE) = c("Intercept", Xnames)
          SE <- as.data.frame(SE)
          df_name <- paste0("group2.row",i)
          se.list2[[df_name]] <- SE

        }
        res <- append(res, se.list2)
      }
    #------------------ SE for group3 ------------------#
      if (!is.null(s2_beta2)) {

      se.list3 <- list()

      for (i in 1:length(s2_beta2)){
        resap.SE.vec <- c()
        for (j in 1:B) {
          mydata.o=mydata
          mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
          mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
          mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
          mydata.o$group=0
          mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1>t0)&(mydata.o$XS2<=t0)]=3
          da3.sap=mydata.o[mydata.o$group==3,]
          G3_model1.sap = loc.fun.ex.gr3.SE(s.seq=s2_beta2, data=da3.sap, t0=t0,  L=L, h=h3.final1,
                                            weight = da3.sap$W.resap)
          G3_coef1.sap = G3_model1.sap$est.mat
          resap.SE.vec = rbind(resap.SE.vec, G3_coef1.sap)
        }
        resap.SE.vec=as.data.frame(resap.SE.vec)
        SE=apply(resap.SE.vec, 2, sd)
        names(SE) = c("Intercept", Xnames)
        SE <- as.data.frame(SE)
        df_name <- paste0("group3.row",i)
        se.list3[[df_name]] <- SE
      }
      res <- append(res, se.list3)
    }
    #------------------ SE for group4 ------------------#
      if (!is.null(s1s2_beta3)) {

      se.list4 <- list()

      for (i in 1:nrow(s1s2_beta3)) {
        resap.SE.vec <- c()
        for (j in 1:B) {
          mydata.o=mydata
          mydata.o$Expwt=rexp(length(mydata.o[,1]), 1)
          mydata.o$what=Wi.FUN(data=mydata.o, t0=t0, tau=L, weight.given=mydata.o$Expwt)
          mydata.o$W.resap=mydata.o$Expwt*mydata.o$what
          mydata.o$group=0
          mydata.o$group[(mydata.o$Y>t0)&(mydata.o$XS1<=t0)&(mydata.o$XS2<=t0)]=4
          da4.sap=mydata.o[mydata.o$group==4,]

          # status4_11=gridSearch(fun=min.BW.cens.ex.gr4.SE,da4=da4.sap,t0=t0,L=L,
          #                       s.seq = s1s2_beta3[i,],
          #                       npar= 2, folds=folds.grid, reps=reps.grid,
          #                       levels = list(a=grid3[[1]]^2, b=grid3[[2]]^2),
          #                       method = gs.method, cl=gs.cl)
          h4.band11_sap=status4_11$minlevels
          h4.final11_sap=(sqrt(h4.band11_sap)/(nrow(da4.sap)^.05))^2
          G4_model11.sap = loc.fun.ex.gr4.SE(s.seq=as.matrix(s1s2_beta3)[1,], data=da4.sap, t0=t0,
                                             L=L, H=rbind(c(h4.final11_sap[1],0), c(0,h4.final11_sap[2])), weight = da4.sap$W.resap)
          G4_coef11.sap = G4_model11.sap$est.mat
          resap.SE.vec = rbind(resap.SE.vec, G4_coef11.sap)
        }
        resap.SE.vec=as.data.frame(resap.SE.vec)
        SE=apply(resap.SE.vec, 2, sd)
        names(SE) = c("Intercept", Xnames)
        SE <- as.data.frame(SE)
        df_name <- paste0("group4.row",i)
        se.list4[[df_name]] <- SE
      }
      res <- append(res, se.list4)
    }
      } else {stop("`SE.gs` must be logical")}

  }
  return(res)
}


