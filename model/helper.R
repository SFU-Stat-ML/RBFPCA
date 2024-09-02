# bisection function
# - recursive implementation of bisection algorithm
bisection <- function(low, high, W, u, phi ){
  mid <- (low+high)/2
  f.low <- rCESS( W, u, low, phi )
  f.mid <- rCESS( W, u, mid, phi )
  f.high <- rCESS( W, u, high, phi )
  
  if( f.low*f.high>0 )
    stop('Invalid endpoint for bisection.')
  
  try({if( low>=high )
    stop('bisection overlap')
    
    if( (abs(f.mid)<1e-10)||((high-low)/2<1e-10) )
      return( mid )
    if( (f.low*f.mid)<0 )
      return( bisection( low, mid, W, u, phi ) )
    if( (f.high*f.mid)<0 )
      return( bisection( mid, high, W, u, phi ) )
  })
  
  stop('bisection flawed')
}

# log-sum-exponential evaluation of log(sum(w))
logsum <- function(logw){
  logmax = max(logw)
  log(sum(exp(logw-logmax)))+logmax
}

# relative conditional effective sample size
rCESS <- function(W, u, a, phi) {
  logw <- a*u          # weight update
  exp(2*logsum(log(W)+logw) - logsum(log(W)+2*logw)) - phi
}

# effective sample size
rESS <- function(logW){
  K <- length(logW)
  logWmax <- max(logW)
  logRESS <- -(2*logWmax + log(sum(exp(2*logW-2*logWmax)))) - log(K)
  return(exp(logRESS))
}

# systematic resampling algorithm
systematic_resample <- function(W){
  K <- length(W)
  U <- runif(1,0,1/K) + 0:(K-1)/K
  W.sum <- cumsum(W)
  N <- rep(NA,K)
  j <- 1
  for( i in 1:K )
  {
    found = F
    while( !found )
    {
      if( U[i]>W.sum[j] )
        j <- j+1
      else
        found = T
    }
    N[i] <- j
  }
  return( N )
}

# remap labels in mixture model to prevent label switching
ReMap.mm <- function(x){
  x.unique <- unique(x)
  x.remap <- vector(mode = "integer", length = length(x))
  for (i in 1:length(x.unique)){
    x.remap[x == x.unique[i]] <- i
  }
  return(x.remap)
}


# The following functions are used to calculate the log marginal likelihood estimates for Gibbs 
BFPCA.robust.MC.beta.fix <- function(data, num.basis, num.PC, n.iter, seed, cov.index, beta){
  # data should be of dimension N.tps * N.subj
  
  set.seed(seed)
  
  # get inputs -----------------------------
  N.subj <- ncol(data) 
  N.tps <- nrow(data) 
  
  J <- num.basis
  K <- num.PC
  
  # Define the hyperparameters -----------------------------
  
  CovFun <- function(s, t, cov.index){
    if (cov.index == 1){
      exp(-3*(t-s)^2)
    } else{
      min(s + 1, t + 1)
    }
  }
  
  ## U_k, L_k ##
  timepts <- seq(-1, 1, length.out = N.tps)
  Sigma.star <- matrix(nrow = N.tps, ncol = N.tps)
  if (cov.index == 0){
    Sigma.star <- cov(t(data))
  } else{
    for (i in 1:N.tps){
      for (j in 1:N.tps){
        Sigma.star[i,j] <- CovFun(s = timepts[i], t = timepts[j], cov.index)
      }
    }
  }
  normalized.p.list <- legendre.polynomials(J-1, normalized=TRUE)
  H <- matrix(nrow = N.tps, ncol = J)
  for (i in 1:length(normalized.p.list)){
    H[,i] <- predict(normalized.p.list[[i]], newdata = timepts)
  }
  HH.inv <- inv(t(H) %*% H)
  xi <- HH.inv %*% t(H) %*% Sigma.star %*% H %*% HH.inv
  Uk <- eigen(xi)$vectors[,1:K]
  Lk <- diag(eigen(xi)$values[1:K])
  UHHU <- t(Uk) %*% t(H) %*% H %*% Uk
  
  ## tau ##
  tau <- 1
  
  ## nu ##
  nu <- 2*K
  
  ## a, b ##
  a <- b <- 10^(-3)
  
  ## theta_0 ##
  theta0 <- rep(0, K)
  
  ##### this part is new -----
  #### start
  ## R ##
  R.diag <- apply(data, 1, function(x) (max(x) - min(x))^2)
  R <- diag(R.diag) 
  
  ## 2r ##
  twor <- N.tps
  
  ## kappa ##
  kappa <- (100/twor) * inv(R) 
  
  ## gamma ##
  gamma <- 10 # each component of vec(D) is given an indep prior N(0, 10^2)
  #### end
  
  # Define the initial values -----------------------------
  
  ## Lambda^-1 ##
  Lambda.inv.path <- list()
  temp1 <- matrix(runif(K^2, min = 0, max = 1)*2-1, ncol=K) 
  Lambda.inv.path[[1]] <- inv(t(temp1) %*% temp1)
  
  ##### this part is new -----
  #### start
  ## Zi ##
  Zi.path <- list()
  Zi.path[[1]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
  for (i in 1:N.subj){
    Zi.path[[1]][i, ] <- mnormt::rmtruncnorm(1, mean = diag(0, N.tps), varcov = diag(1, N.tps), 
                                             lower = rep(0, N.tps), upper = rep(Inf, N.tps))
  }
  
  ## vec(D) ##
  D.path <- list()
  vecD.diag <- rnorm(n = N.tps, mean = 0, sd = gamma)
  D.path[[1]] <- vecD.diag 
  
  ## Sigma^-1 ##
  Sigma.inv.path <- list()
  temp2 <- matrix(runif(N.tps^2, min = 0, max = 1)*2-1, ncol=N.tps)
  Sigma.inv.path[[1]] <- inv(t(temp2) %*% temp2)
  
  #### end
  
  # Run Gibbs sampler -----------------------------
  
  tic()
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  for (iter in 2:n.iter){
    
    Lambda.inv.path[[iter]] <- matrix(NA, nrow = K, ncol = K)
    Zi.path[[iter]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
    
    ## sample Lambda^-1 from a Wishart distribution
    Lambda.inv.df <- nu + N.subj + 1
    Lambda.inv.p2 <- matrix(0, nrow = K, ncol = K)
    for (i in 1:N.subj){
      Lambda.inv.p2 <- Lambda.inv.p2 + 
        (beta[i,]) %*% t(beta[i,])
    }
    Lambda.inv.mat <- inv(Lk + Lambda.inv.p2)
    Lambda.inv.path[[iter]] <- MCMCpack::rwish(v = Lambda.inv.df, S = Lambda.inv.mat)
    
    ##### this part is new -----
    #### start
    ## sample Zi from a T-dim Gaussian distribution
    Ai.temp <- diag(D.path[[iter-1]]) %*% Sigma.inv.path[[iter-1]]
    Ai <- diag(N.tps) + Ai.temp %*% diag(D.path[[iter-1]])
    Zi.path[[iter]] <- foreach(i=1:N.subj, .combine='rbind', .multicombine=TRUE) %dopar% {
      ai <- Ai.temp %*% (data[, i] - H %*% Uk %*% beta[i,])
      Zi.mean <- inv(Ai) %*% ai
      Zi.varcov <- inv(Ai)
      mnormt::rmtruncnorm(1, mean = Zi.mean, varcov = Zi.varcov,
                          lower = rep(0, N.tps), upper = rep(Inf, N.tps))
    }
    
    ## sample vec(D) from a T-dim Gaussian distribution
    add <- function(x) Reduce("+", x)
    comb <- function(...) {
      mapply('rbind', ..., SIMPLIFY=FALSE)
    }
    res <- foreach(i=1:N.subj, .combine='comb', .multicombine=TRUE) %dopar% {
      B.temp <- diag(Zi.path[[iter]][i, ]) %*% Sigma.inv.path[[iter-1]] %*% diag(Zi.path[[iter]][i, ])
      b.temp <- diag(Zi.path[[iter]][i, ]) %*% Sigma.inv.path[[iter-1]] %*%
        (data[, i] - H %*% Uk %*% beta[i,]) %>% as.numeric()
      list(B.temp = B.temp, b.temp = b.temp)
    }
    number_of_chunks = nrow(res$B.temp) / N.tps
    B.temp.list = lapply(split(res$B.temp, rep(1:number_of_chunks, each = NROW(res$B.temp)/number_of_chunks)),
                         function(a) matrix(a, ncol = NCOL(res$B.temp)))
    B.temp <- add(B.temp.list) + diag(x = 1/gamma, nrow = N.tps)
    b.temp <- colSums(res$b.temp)
    vecD.mean <- inv(B.temp) %*% b.temp
    vecD.varcov <- inv(B.temp)
    vecD.diag <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
    D.path[[iter]] <- vecD.diag
    
    ## sample Sigma^(-1) from a Wishart distribution
    Sigma.inv.df <- twor + N.subj
    Sigma.inv.mat <- inv(2*kappa)
    for (i in 1:N.subj){
      temp <- data[, i] - H %*% Uk %*% beta[i,] - diag(D.path[[iter]]) %*% Zi.path[[iter]][i, ]
      Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
    }
    Sigma.inv.path[[iter]] <- MCMCpack::rwish(v = Sigma.inv.df, S = inv(Sigma.inv.mat))
    
    #### end
    
    print(paste0("iteration ", iter, " is done"))
    
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  toc()
  
  results <- list(Lambda.inv = Lambda.inv.path,            
                  Zi = Zi.path, 
                  D = D.path,
                  Sigma.inv = Sigma.inv.path,
                  H = H, Uk = Uk, Lk = Lk)  
  
  return(results)
  
}

BFPCA.robust.MC.beta.D.fix <- function(data, num.basis, num.PC, n.iter, seed, cov.index, beta, D){
  # data should be of dimension N.tps * N.subj
  
  set.seed(seed)
  
  # get inputs -----------------------------
  N.subj <- ncol(data) 
  N.tps <- nrow(data) 
  
  J <- num.basis
  K <- num.PC
  
  # Define the hyperparameters -----------------------------
  
  CovFun <- function(s, t, cov.index){
    if (cov.index == 1){
      exp(-3*(t-s)^2)
    } else{
      min(s + 1, t + 1)
    }
  }
  
  ## U_k, L_k ##
  timepts <- seq(-1, 1, length.out = N.tps)
  Sigma.star <- matrix(nrow = N.tps, ncol = N.tps)
  if (cov.index == 0){
    Sigma.star <- cov(t(data))
  } else{
    for (i in 1:N.tps){
      for (j in 1:N.tps){
        Sigma.star[i,j] <- CovFun(s = timepts[i], t = timepts[j], cov.index)
      }
    }
  }
  normalized.p.list <- legendre.polynomials(J-1, normalized=TRUE)
  H <- matrix(nrow = N.tps, ncol = J)
  for (i in 1:length(normalized.p.list)){
    H[,i] <- predict(normalized.p.list[[i]], newdata = timepts)
  }
  HH.inv <- inv(t(H) %*% H)
  xi <- HH.inv %*% t(H) %*% Sigma.star %*% H %*% HH.inv
  Uk <- eigen(xi)$vectors[,1:K]
  Lk <- diag(eigen(xi)$values[1:K])
  UHHU <- t(Uk) %*% t(H) %*% H %*% Uk
  
  ## tau ##
  tau <- 1
  
  ## nu ##
  nu <- 2*K
  
  ## a, b ##
  a <- b <- 10^(-3)
  
  ## theta_0 ##
  theta0 <- rep(0, K)
  
  ##### this part is new -----
  #### start
  ## R ##
  R.diag <- apply(data, 1, function(x) (max(x) - min(x))^2)
  R <- diag(R.diag) 
  
  ## 2r ##
  twor <- N.tps
  
  ## kappa ##
  kappa <- (100/twor) * inv(R) 
  
  ## gamma ##
  gamma <- 10 # each component of vec(D) is given an indep prior N(0, 10^2)
  #### end
  
  # Define the initial values -----------------------------
  
  ## Lambda^-1 ##
  Lambda.inv.path <- list()
  temp1 <- matrix(runif(K^2, min = 0, max = 1)*2-1, ncol=K) 
  Lambda.inv.path[[1]] <- inv(t(temp1) %*% temp1)
  
  ##### this part is new -----
  #### start
  ## Zi ##
  Zi.path <- list()
  Zi.path[[1]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
  for (i in 1:N.subj){
    Zi.path[[1]][i, ] <- mnormt::rmtruncnorm(1, mean = diag(0, N.tps), varcov = diag(1, N.tps), 
                                             lower = rep(0, N.tps), upper = rep(Inf, N.tps))
  }
  
  ## Sigma^-1 ##
  Sigma.inv.path <- list()
  temp2 <- matrix(runif(N.tps^2, min = 0, max = 1)*2-1, ncol=N.tps)
  Sigma.inv.path[[1]] <- inv(t(temp2) %*% temp2)
  
  #### end
  
  # Run Gibbs sampler -----------------------------
  
  tic()
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  for (iter in 2:n.iter){
    
    Lambda.inv.path[[iter]] <- matrix(NA, nrow = K, ncol = K)
    Zi.path[[iter]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
    
    ## sample Lambda^-1 from a Wishart distribution
    Lambda.inv.df <- nu + N.subj + 1
    Lambda.inv.p2 <- matrix(0, nrow = K, ncol = K)
    for (i in 1:N.subj){
      Lambda.inv.p2 <- Lambda.inv.p2 + 
        (beta[i,]) %*% t(beta[i,])
    }
    Lambda.inv.mat <- inv(Lk + Lambda.inv.p2)
    Lambda.inv.path[[iter]] <- MCMCpack::rwish(v = Lambda.inv.df, S = Lambda.inv.mat)
    
    ##### this part is new -----
    #### start
    ## sample Zi from a T-dim Gaussian distribution
    Ai.temp <- diag(D) %*% Sigma.inv.path[[iter-1]]
    Ai <- diag(N.tps) + Ai.temp %*% diag(D)
    Zi.path[[iter]] <- foreach(i=1:N.subj, .combine='rbind', .multicombine=TRUE) %dopar% {
      ai <- Ai.temp %*% (data[, i] - H %*% Uk %*% beta[i,])
      Zi.mean <- inv(Ai) %*% ai
      Zi.varcov <- inv(Ai)
      mnormt::rmtruncnorm(1, mean = Zi.mean, varcov = Zi.varcov,
                          lower = rep(0, N.tps), upper = rep(Inf, N.tps))
    }
    
    ## sample Sigma^(-1) from a Wishart distribution
    Sigma.inv.df <- twor + N.subj
    Sigma.inv.mat <- inv(2*kappa)
    for (i in 1:N.subj){
      temp <- data[, i] - H %*% Uk %*% beta[i,] - diag(D) %*% Zi.path[[iter]][i, ]
      Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
    }
    Sigma.inv.path[[iter]] <- MCMCpack::rwish(v = Sigma.inv.df, S = inv(Sigma.inv.mat))
    
    #### end
    
    print(paste0("iteration ", iter, " is done"))
    
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  toc()
  
  results <- list(Lambda.inv = Lambda.inv.path,            
                  Zi = Zi.path, 
                  Sigma.inv = Sigma.inv.path,
                  H = H, Uk = Uk, Lk = Lk)  
  
  return(results)
  
}

BFPCA.robust.MC.beta.D.Sigma.inv.fix <- function(data, num.basis, num.PC, n.iter, seed, cov.index, beta, D, Sigma.inv){
  # data should be of dimension N.tps * N.subj
  
  set.seed(seed)
  
  # get inputs -----------------------------
  N.subj <- ncol(data) 
  N.tps <- nrow(data) 
  
  # Posterior computation is done independently for each pair (J,K), 1 ≤ K ≤ J ≤ Jmax.
  J <- num.basis
  K <- num.PC
  
  
  # Define the hyperparameters -----------------------------
  
  CovFun <- function(s, t, cov.index){
    if (cov.index == 1){
      exp(-3*(t-s)^2)
    } else{
      min(s + 1, t + 1)
    }
  }
  
  ## U_k, L_k ##
  timepts <- seq(-1, 1, length.out = N.tps)
  Sigma.star <- matrix(nrow = N.tps, ncol = N.tps)
  if (cov.index == 0){
    Sigma.star <- cov(t(data))
  } else{
    for (i in 1:N.tps){
      for (j in 1:N.tps){
        Sigma.star[i,j] <- CovFun(s = timepts[i], t = timepts[j], cov.index)
      }
    }
  }
  normalized.p.list <- legendre.polynomials(J-1, normalized=TRUE)
  H <- matrix(nrow = N.tps, ncol = J)
  for (i in 1:length(normalized.p.list)){
    H[,i] <- predict(normalized.p.list[[i]], newdata = timepts)
  }
  HH.inv <- inv(t(H) %*% H)
  xi <- HH.inv %*% t(H) %*% Sigma.star %*% H %*% HH.inv
  Uk <- eigen(xi)$vectors[,1:K]
  Lk <- diag(eigen(xi)$values[1:K])
  UHHU <- t(Uk) %*% t(H) %*% H %*% Uk
  
  ## tau ##
  tau <- 1
  
  ## nu ##
  nu <- 2*K
  
  ## a, b ##
  a <- b <- 10^(-3)
  
  ## theta_0 ##
  theta0 <- rep(0, K)
  
  ##### this part is new -----
  #### start
  ## R ##
  R.diag <- apply(data, 1, function(x) (max(x) - min(x))^2)
  R <- diag(R.diag) 
  
  ## 2r ##
  twor <- N.tps
  
  ## kappa ##
  kappa <- (100/twor) * inv(R) 
  
  ## gamma ##
  gamma <- 10 # each component of vec(D) is given an indep prior N(0, 10^2)
  #### end
  
  # Define the initial values -----------------------------
  
  ## Lambda^-1 ##
  Lambda.inv.path <- list()
  temp1 <- matrix(runif(K^2, min = 0, max = 1)*2-1, ncol=K) 
  Lambda.inv.path[[1]] <- inv(t(temp1) %*% temp1)
  
  ##### this part is new -----
  #### start
  ## Zi ##
  Zi.path <- list()
  Zi.path[[1]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
  for (i in 1:N.subj){
    Zi.path[[1]][i, ] <- mnormt::rmtruncnorm(1, mean = diag(0, N.tps), varcov = diag(1, N.tps), 
                                             lower = rep(0, N.tps), upper = rep(Inf, N.tps))
  }
  
  #### end
  
  # Run Gibbs sampler -----------------------------
  
  tic()
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  for (iter in 2:n.iter){
    
    Lambda.inv.path[[iter]] <- matrix(NA, nrow = K, ncol = K)
    Zi.path[[iter]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
    
    ## sample Lambda^-1 from a Wishart distribution
    Lambda.inv.df <- nu + N.subj + 1
    Lambda.inv.p2 <- matrix(0, nrow = K, ncol = K)
    for (i in 1:N.subj){
      Lambda.inv.p2 <- Lambda.inv.p2 + 
        (beta[i,]) %*% t(beta[i,])
    }
    Lambda.inv.mat <- inv(Lk + Lambda.inv.p2)
    Lambda.inv.path[[iter]] <- MCMCpack::rwish(v = Lambda.inv.df, S = Lambda.inv.mat)
    
    ##### this part is new -----
    #### start
    ## sample Zi from a T-dim Gaussian distribution
    Ai.temp <- diag(D) %*% Sigma.inv
    Ai <- diag(N.tps) + Ai.temp %*% diag(D)
    Zi.path[[iter]] <- foreach(i=1:N.subj, .combine='rbind', .multicombine=TRUE) %dopar% {
      ai <- Ai.temp %*% (data[, i] - H %*% Uk %*% beta[i,])
      Zi.mean <- inv(Ai) %*% ai
      Zi.varcov <- inv(Ai)
      mnormt::rmtruncnorm(1, mean = Zi.mean, varcov = Zi.varcov,
                          lower = rep(0, N.tps), upper = rep(Inf, N.tps))
    }
    
    #### end
    
    print(paste0("iteration ", iter, " is done"))
    
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  toc()
  
  results <- list(Lambda.inv = Lambda.inv.path,            
                  Zi = Zi.path, 
                  H = H, Uk = Uk, Lk = Lk)  
  
  return(results)
  
}


checkTmvArgs <- function(mean, sigma, lower, upper){
  if (is.null(lower) || any(is.na(lower))) 
    stop(sQuote("lower"), " not specified or contains NA")
  if (is.null(upper) || any(is.na(upper))) 
    stop(sQuote("upper"), " not specified or contains NA")
  if (!is.numeric(mean) || !is.vector(mean)) 
    stop(sQuote("mean"), " is not a numeric vector")
  if (is.null(sigma) || any(is.na(sigma))) 
    stop(sQuote("sigma"), " not specified or contains NA")
  
  if (!is.matrix(sigma)) {
    sigma <- as.matrix(sigma)
  }
  
  if (NCOL(lower) != NCOL(upper)) {
    stop("lower and upper have non-conforming size")
  }
  
  checkSymmetricPositiveDefinite(sigma)
  
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  
  if (length(lower) != length(mean) || length(upper) != length(mean)) {
    stop("mean, lower and upper must have the same length")
  }
  
  if (any(lower>=upper)) {
    stop("lower bound should be strictly less than the upper bound (lower<upper)")
  }
  
  # checked arguments
  cargs <- list(mean=mean, sigma=sigma, lower=lower, upper=upper)
  return(cargs)
}

checkSymmetricPositiveDefinite <- function(x, name="sigma") {
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a symmetric matrix", name))
  }
  
  if (NROW(x) != NCOL(x)) {
    stop(sprintf("%s must be a square matrix", name))
  }
  
  if (any(diag(x) <= 0)) {
    stop(sprintf("%s all diagonal elements must be positive", name))
  }
  
  if (det(x) <= 0) {
    stop(sprintf("%s must be positive definite", name))
  }
}

mydtmvnorm <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                       lower = rep( -Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), 
                       log = FALSE, margin=NULL){
  # check of standard tmvnorm arguments
  cargs <- checkTmvArgs(mean=mean, sigma=sigma, lower=lower, upper=upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  # Check of optional argument "margin"
  if (!is.null(margin)) {
    if (!length(margin) %in% c(1, 2))
      stop("Length of margin must be either 1 (one-dimensional marginal density) or 2 (bivariate marginal density).")
    if (any(margin <= 0) || any(margin > length(mean))) {
      stop("All elements in margin must be in 1..length(mean).")	
    }
    # one-dimensional marginal density f_{n}(x_n)
    if (length(margin) == 1) {
      return(dtmvnorm.marginal(xn=x, n=margin, mean = mean, sigma = sigma, lower = lower, upper = upper, log = log))		
    }
    # for bivariate marginal density f_{q,r}(x_q, x_r) we need q <> r and "x" as (n x 2) matrix
    if (length(margin) == 2) {
      if(margin[1] == margin[2])	
        stop("Two different margins needed for bivariate marginal density.")
      if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
      }  
      if(!is.matrix(x) || ncol(x) != 2)
        stop("For bivariate marginal density x must be either a (n x 2) matrix or a vector of length 2.")  
      # bivariate marginal density f_{q,r}(x_q, x_r)	
      return(dtmvnorm.marginal2(xq=x[,1], xr=x[,2], q=margin[1], r=margin[2], mean = mean, sigma = sigma, lower = lower, upper = upper, log = log))	  
    }	
  }
  
  # Check of range of x
  if (min(x) < min(lower)){
    x <- abs(x)
  }
  
  # Check of additional inputs like x
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  
  # Anzahl der Beobachtungen
  T <- nrow(x)
  
  # check for each row if in support region
  insidesupportregion <- logical(T)
  for (i in 1:T)
  {
    insidesupportregion[i] = all(x[i,] >= lower & x[i,] <= upper & !any(is.infinite(x)))
  }
  
  if(log) {
    # density value for points inside the support region
    dvin <- dmvnorm(x, mean=mean, sigma=sigma, log=TRUE) - log(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)) 
    # density value for points outside the support region
    dvout <- -Inf
  } else {
    dvin <- dmvnorm(x, mean=mean, sigma=sigma, log=FALSE) / pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
    dvout <- 0
  }
  
  
  f <- ifelse(insidesupportregion, dvin, dvout)
  return(f)
}

mydwish <- function(x, nu, S){
  # give density in normal scale
  det(S)^(nu/2) * det(x)^((nu-nrow(x)-1)/2) * exp(-tr(S %*% x)/2)
}

mydwish2 <- function(x, nu, S){
  # give density ini log scale
  log(det(S)^((nu/2)/5)) * 5 + log(det(x)^((nu-nrow(x)-1)/2/4)) * 4 - tr(S %*% x)/2
}

logL.fun <- function(data, beta, Lambda.inv, Zi, D, Sigma.inv){
  
  logL <- 0
  for (j in 1:ncol(data)){
    mvnorm.mean <- BFPCA.robust.res$H %*% BFPCA.robust.res$Uk %*% beta[j, ] + diag(D) %*% Zi[j, ]
    result <- mvtnorm::dmvnorm(x = data[, j], mean = mvnorm.mean, sigma = inv(Sigma.inv), log = TRUE)
    logL <- logL + result
  }
  return(logL)
}

log.prior.fun <- function(data, beta, Lambda.inv, Zi, D, Sigma.inv){
  
  log.prior <- 0
  for (j in 1:ncol(data)){
    temp <- mvtnorm::dmvnorm(beta[j,], mean = rep(0, ncol(beta)), sigma = inv(Lambda.inv), log = TRUE) + 
      mydtmvnorm(Zi[j,], mean = rep(0, ncol(Zi)), sigma = diag(x = 1, ncol = ncol(Zi), nrow = ncol(Zi)),
                 lower = rep(0, ncol(Zi)), upper = rep(Inf, ncol(Zi)), log = TRUE)
    log.prior <- log.prior + temp
  }
  
  R.diag <- apply(data, 1, function(x) (max(x) - min(x))^2)
  R <- diag(R.diag) 
  twor <- K.FPCA
  kappa <- ((100/twor)*inv(R))[1:K.FPCA, 1:K.FPCA]
  
  log.prior.all <- log.prior + log(MCMCpack::dwish(Lambda.inv, v = 2*nrow(Lambda.inv), S = inv(BFPCA.robust.res$Lk))) +
    sum(sapply(D, dnorm, mean = 0, sd = gamma, log = TRUE)) +
    log(mydwish(x = Sigma.inv[1:K.FPCA, 1:K.FPCA], nu = twor, S = 2*kappa))
  
  return(log.prior.all)
}


pi.beta <- function(data, beta, Lambda.inv, Sigma.inv, Uk, H){
  result <- 0
  for (i in 1:nrow(beta)){
    beta.mean1 <- inv(t(Uk) %*% t(H) %*% Sigma.inv %*% H %*% Uk + Lambda.inv) 
    beta.mean2 <- t(Uk) %*% t(H) %*% Sigma.inv %*% data[,i] 
    beta.mean <- beta.mean1 %*% beta.mean2
    beta.variance <- beta.mean1
    result <- result + mvtnorm::dmvnorm(x = beta[i,], mean = beta.mean, sigma = beta.variance, log = FALSE)
  }
  return(result)
}

pi.D <- function(data, beta, Lambda.inv, Zi, Sigma.inv, D, H, Uk){
  N.subj <- ncol(data)
  N.tps <- nrow(data)
  add <- function(x) Reduce("+", x)
  comb <- function(...) {
    mapply('rbind', ..., SIMPLIFY=FALSE)
  }
  res <- foreach(i=1:N.subj, .combine='comb', .multicombine=TRUE) %dopar% {
    B.temp <- diag(Zi[i, ]) %*% Sigma.inv %*% diag(Zi[i, ])
    b.temp <- diag(Zi[i, ]) %*% Sigma.inv %*%
      (data[, i] - H %*% Uk %*% beta[i,]) %>% as.numeric()
    list(B.temp = B.temp, b.temp = b.temp)
  }
  number_of_chunks = nrow(res$B.temp) / N.tps
  B.temp.list = lapply(split(res$B.temp, rep(1:number_of_chunks, each = NROW(res$B.temp)/number_of_chunks)),
                       function(a) matrix(a, ncol = NCOL(res$B.temp)))
  B.temp <- add(B.temp.list) + diag(x = 1/gamma, nrow = N.tps)
  b.temp <- colSums(res$b.temp)
  vecD.mean <- inv(B.temp) %*% b.temp
  vecD.varcov <- inv(B.temp)
  result <- mnormt::dmnorm(x = as.numeric(D), mean = as.numeric(vecD.mean), varcov = vecD.varcov, log = TRUE)
  print(result)
  return(result)
}

pi.Zi <- function(data, beta, Lambda.inv, Zi, Sigma.inv, D, H, Uk){
  N.tps <- nrow(data)
  Ai.temp <- diag(D) %*% Sigma.inv
  Ai <- diag(N.tps) + Ai.temp %*% diag(D)
  result <- 0
  for (i in 1:nrow(beta)){
    ai <- Ai.temp %*% (data[, i] - H %*% Uk %*% beta[i,])
    Zi.mean <- inv(Ai) %*% ai
    Zi.varcov <- inv(Ai)
    result <- result + mvtnorm::dmvnorm(Zi[i,], mean = as.numeric(Zi.mean), sigma = Zi.varcov, log = TRUE)
  }
  return(result)
}

pi.Sigma.inv <- function(data, beta, Zi, Sigma.inv, D, H, Uk){
  N.subj <- ncol(data)
  R.diag <- apply(data, 1, function(x) (max(x) - min(x))^2)
  R <- diag(R.diag)
  twor <- nrow(data)
  kappa <- ((100/twor)*inv(R))
  Sigma.inv.df <- twor + N.subj
  Sigma.inv.mat <- inv(2*kappa)
  for (i in 1:N.subj){
    temp <- data[, i] - H %*% Uk %*% beta[i,] - diag(D) %*% Zi[i, ]
    Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
  }
  result <- mydwish2(x = Sigma.inv, nu = Sigma.inv.df, S = inv(Sigma.inv.mat))
  print(result)
  return(result)
}

pi.Lambda.inv <- function(data, beta, Lambda.inv, Sigma.inv, Lk){
  K <- nrow(Lambda.inv)
  nu <- 2*K
  N.subj <- nrow(beta)
  Lambda.inv.df <- nu + N.subj + 1
  Lambda.inv.p2 <- matrix(0, nrow = K, ncol = K)
  for (i in 1:N.subj){
    Lambda.inv.p2 <- Lambda.inv.p2 +
      (beta[i,]) %*% t(beta[i,])
  }
  Lambda.inv.mat <- inv(Lk + Lambda.inv.p2)
  result <- mydwish2(x = Lambda.inv, nu = Lambda.inv.df, S = Lambda.inv.mat)
  print(result)
  return(result)
}


