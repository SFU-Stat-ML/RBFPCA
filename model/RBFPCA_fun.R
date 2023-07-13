RBFPCA <- function(data, num.basis, num.PC, n.iter, seed, cov.index){
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
  
  ## beta_{i,k} ##
  beta.path <- list()
  beta.path[[1]] <- matrix(0, nrow = N.subj, ncol = K)
  
  ## Lambda^-1 ##
  Lambda.inv.path <- list()
  temp1 <- matrix(runif(K^2, min = 0, max = 1)*2-1, ncol=K) 
  Lambda.inv.path[[1]] <- inv(t(temp1) %*% temp1)
  
  ##### this part is for skew elliptical distribution -----
  #### start
  ## Zi ##
  Zi.path <- list()
  Zi.path[[1]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
  for (i in 1:N.subj){
    Zi.path[[1]][i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps,
                                          mean = rep(0, N.tps), sigma = diag(1, N.tps),
                                          lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
  }
  
  ## vec(D) ##
  D.path <- list()
  vecD.diag <- rnorm(n = N.tps, mean = 0, sd = gamma)
  D.path[[1]] <- vecD.diag 
  
  ## Sigma^-1 ##
  Sigma.inv.path <- list()
  temp2 <- matrix(runif(N.tps^2, min = 0, max = 10), ncol=N.tps)
  Sigma.inv.path[[1]] <- inv(t(temp2) %*% temp2)
  

  # Run Gibbs sampler -----------------------------
  
  tic()
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  for (iter in 2:n.iter){
    beta.path[[iter]] <- matrix(NA, nrow = N.subj, ncol = K)
    Lambda.inv.path[[iter]] <- matrix(NA, nrow = K, ncol = K)
    Zi.path[[iter]] <- matrix(data = NA, nrow = N.subj, ncol = N.tps)

    ## for each subject, sample beta_{i,k} from a K-dim Gaussian distribution
    beta.path[[iter]] <- foreach(i=1:N.subj, .combine='rbind', .multicombine=TRUE) %dopar% {
      beta.mean1 <- inv(t(Uk) %*% t(H) %*% Sigma.inv.path[[iter-1]] %*% H %*% Uk + Lambda.inv.path[[iter-1]]) 
      beta.mean2 <- t(Uk) %*% t(H) %*% Sigma.inv.path[[iter-1]] %*% data[,i] 
      beta.mean <- beta.mean1 %*% beta.mean2
      beta.variance <- beta.mean1
      MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
    }

    ## sample Lambda^-1 from a Wishart distribution
    Lambda.inv.df <- nu + N.subj + 1
    Lambda.inv.p2 <- matrix(0, nrow = K, ncol = K)
    for (i in 1:N.subj){
      Lambda.inv.p2 <- Lambda.inv.p2 + 
        (beta.path[[iter]][i,]) %*% t(beta.path[[iter]][i,])
    }
    Lambda.inv.mat <- inv(Lk + Lambda.inv.p2)
    Lambda.inv.path[[iter]] <- MCMCpack::rwish(v = Lambda.inv.df, S = Lambda.inv.mat)
    
    ## sample Zi from a T-dim Gaussian distribution
    Ai.temp <- diag(D.path[[iter-1]]) %*% Sigma.inv.path[[iter-1]]
    Ai <- diag(N.tps) + Ai.temp %*% diag(D.path[[iter-1]])
    Zi.path[[iter]] <- foreach(i=1:N.subj, .combine='rbind', .multicombine=TRUE) %dopar% {
      ai <- Ai.temp %*% (data[, i] - H %*% Uk %*% beta.path[[iter]][i,])
      Zi.mean <- inv(Ai) %*% ai
      Zi.varcov <- inv(Ai)
      tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov),
                       mean = as.numeric(Zi.mean), sigma = Zi.varcov,
                       lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
    }
    
    ## sample vec(D) from a T-dim Gaussian distribution
    add <- function(x) Reduce("+", x)
    comb <- function(...) {
      mapply('rbind', ..., SIMPLIFY=FALSE)
    }
    res <- foreach(i=1:N.subj, .combine='comb', .multicombine=TRUE) %dopar% {
      B.temp <- diag(Zi.path[[iter]][i, ]) %*% Sigma.inv.path[[iter-1]] %*% diag(Zi.path[[iter]][i, ])
      b.temp <- diag(Zi.path[[iter]][i, ]) %*% Sigma.inv.path[[iter-1]] %*%
        (data[, i] - H %*% Uk %*% beta.path[[iter]][i,]) %>% as.numeric()
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
      temp <- data[, i] - H %*% Uk %*% beta.path[[iter]][i,] - diag(D.path[[iter]]) %*% Zi.path[[iter]][i, ]
      Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
    }
    Sigma.inv.path[[iter]] <- MCMCpack::rwish(v = Sigma.inv.df, S = inv(Sigma.inv.mat))
    
    print(paste0("iteration ", iter, " is done"))
    
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  toc()
  
  results <- list(beta = beta.path,          
                  Lambda.inv = Lambda.inv.path,            
                  Zi = Zi.path, 
                  D = D.path,
                  Sigma.inv = Sigma.inv.path,
                  H = H, Uk = Uk, Lk = Lk)  
  
  return(results)
  
}

check.traceplot.RBFPCA <- function(RBFPCA.result, iter.burnin, n.burnin = 1, n.iter, 
                                   i.plot, t.plot, k1.plot, k2.plot){
  
  # get inputs
  beta.path <- RBFPCA.result$beta
  Lambda.inv.path <- RBFPCA.result$Lambda.inv
  Zi.path <- RBFPCA.result$Zi
  D.path <- RBFPCA.result$D
  Sigma.inv.path <- RBFPCA.result$Sigma.inv
  
  # Get trace plots -----------------------------
  
  par(mfrow = c(2, 3))
  
  ## for beta_{i,k}
  beta.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    beta.df[iter] <- beta.path[[iter]][i.plot,k1.plot]
  }
  plot(beta.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "beta",
       main = paste0("Evaluated at i = ", i.plot, " and k = ", k1.plot))
  abline(v = iter.burnin, col = "blue")
  
  ## for Lambda^-1
  Lambda.inv.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    Lambda.inv.df[iter] <- Lambda.inv.path[[iter]][k1.plot, k2.plot]
  }
  plot(Lambda.inv.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "Lambda.inv",
       main = paste0("Evaluated at k1 = ", k1.plot, " and k2 = ", k2.plot))
  abline(v = iter.burnin, col = "blue")
  
  ## for Zi
  Zi.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    Zi.df[iter] <- Zi.path[[iter]][i.plot, t.plot]
  }
  plot(Zi.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "Zi",
       main = paste0("Evaluated at i = ", i.plot, " and t = ", t.plot))
  abline(v = iter.burnin, col = "blue")
  
  ## for D
  D.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    D.df[iter] <- D.path[[iter]][t.plot]
  }
  plot(D.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "D",
       main = paste0("Evaluated at t = ", t.plot))
  abline(v = iter.burnin, col = "blue")
  
  ## for Sigma^-1
  Sigma.inv.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    Sigma.inv.df[iter] <- Sigma.inv.path[[iter]][t.plot, t.plot]
  }
  plot(Sigma.inv.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "Sigma.inv",
       main = paste0("Evaluated at t = ", t.plot))
  abline(v = iter.burnin, col = "blue")
  
  par(mfrow = c(1, 1))
  
}


RBFPCA.sparse <- function(data.table, data, num.basis, num.PC, n.iter, seed, cov.index, search.bw){
  
  set.seed(seed)
  
  N.subj <- length(data$x.demean) 
  J <- num.basis
  K <- num.PC
  
  # Define the hyperparameters -----------------------------
  
  CovFun <- function(s, t, cov.index){
    if (cov.index == 1){
      exp(-3*(t-s)^2)
    } else if (cov.index == 2){
      min(s + 1, t + 1)
    } else{
      exp(-1*(t-s)^2)
    }
  }
  
  # Define the initial values for parameters whose sizes don't depend on n_i -----------------------------
  
  ## beta_{i,k} ##
  beta.path <- list()
  beta.path[[1]] <- matrix(0, nrow = N.subj, ncol = K)
  
  ## Lambda^-1 ##
  Lambda.inv.path <- list()
  temp1 <- matrix(runif(K^2, min = 0, max = 1)*2-1, ncol=K) 
  Lambda.inv.path[[1]] <- inv(t(temp1) %*% temp1)
  
  # Initialize the containers for parameters whose sizes depend on n_i
  Zi.path <- list()
  Zi.path[[1]] <- list()
  D.path <- list()
  D.path[[1]] <- list()
  Sigma.inv.path <- list()
  Sigma.inv.path[[1]] <- list()
  
  H.list <- list()
  Uk.list <- list()
  Lk.list <- list()
  kappa.list <- list()
  
  get.corner <- function(x){
    c.end <- ncol(x)
    r.end <- nrow(x)
    a = x[1,1]
    b = x[1,c.end]
    c = x[r.end,1]
    d = x[r.end,c.end]
    matrix(c(a,b,c,d), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  for (i in 1:N.subj){
    
    N.tps <- length(data$pp[[i]]) # number of observations for this subject
    this.id <- names(data$pp)[[i]]  # id for this subject
    
    ## U_k, L_k ##
    timepts <- seq(-1, 1, length.out = N.tps)
    Sigma.star <- matrix(nrow = N.tps, ncol = N.tps)
    if (cov.index == 0){
      if (N.tps == 2){
        PACE.res <- fdapace::FPCA(Ly = data$x, Lt = data$pp, list(FVEthreshold = 0.99, nRegGrid = N.tps+1))
        Sigma.star <- get.corner(PACE.res$fittedCov)
      } else{
        PACE.res <- fdapace::FPCA(Ly = data$x, Lt = data$pp, list(FVEthreshold = 0.99, nRegGrid = N.tps))
        Sigma.star <- PACE.res$fittedCov
      }
    } else{
      for (ii in 1:N.tps){
        for (jj in 1:N.tps){
          Sigma.star[ii,jj] <- CovFun(s = timepts[ii], t = timepts[jj], cov.index)
        }
      }
    }
    normalized.p.list <- legendre.polynomials(J-1, normalized=TRUE)
    H <- matrix(nrow = N.tps, ncol = J)
    for (ii in 1:length(normalized.p.list)){
      H[,ii] <- predict(normalized.p.list[[ii]], newdata = timepts)
    }
    HH.inv <- inv(t(H) %*% H)
    xi <- HH.inv %*% t(H) %*% Sigma.star %*% H %*% HH.inv
    Uk <- eigen(xi)$vectors[,1:K]
    Lk <- diag(eigen(xi)$values[1:K])
    UHHU <- t(Uk) %*% t(H) %*% H %*% Uk
    
    H.list[[i]] <- H
    Uk.list[[i]] <- Uk
    Lk.list[[i]] <- Lk
    
    ## tau ##
    tau <- 1
    
    ## nu ##
    nu <- 2*K
    
    ## a, b ##
    a = b <- 10^(-3)
    
    ## theta_0 ##
    theta0 <- rep(0, K)
    
    ## R ##
    R.diag.fun <- function(x){
      (max(x) - min(x))^2
    }
    R.diag <- c()
    for (tp in 1:N.tps){
      target.pp <- data$pp[[i]][tp]
      search.data <- data.table %>%
        filter(id != this.id) %>%
        mutate(diff = abs(pp-target.pp)) %>%
        arrange(diff) %>% 
        slice_head(n = search.bw)
      R.diag[tp] <- R.diag.fun(search.data$x.demean)
    }
    R <- diag(R.diag)
    
    ## 2r ##
    twor <- N.tps
    
    ## kappa ##
    kappa <- 100 / twor * inv(R)
    kappa.list[[i]] <- kappa
    
    ## gamma ##
    gamma <- 10 # each component of vec(D) is given an indep prior N(0, 10^2)
    #### end
    
    # Define the initial values for parameters whose sizes depend on n_i -----------------------------
    
    ## Zi ##
    Zi.path[[1]][[this.id]] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps,
                                                mean = rep(0, N.tps), sigma = diag(1, N.tps),
                                                lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
    
    ## vec(D) ##
    vecD.diag <- rnorm(n = N.tps, mean = 0, sd = gamma)
    D.path[[1]][[this.id]] <- vecD.diag 
    
    ## Sigma^-1 ##
    temp2 <- matrix(runif(N.tps^2, min = 0, max = 1)*2-1, ncol=N.tps)
    Sigma.inv.path[[1]][[this.id]] <- inv(t(temp2) %*% temp2)
    
    print(paste0("preprocessing: ", i, "/", N.subj, " is done"))
  }
  
  # Run Gibbs sampler -----------------------------
  
  tic()
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  for (iter in 2:n.iter){
    
    beta.path[[iter]] <- matrix(NA, nrow = N.subj, ncol = K)
    Lambda.inv.path[[iter]] <- matrix(NA, nrow = K, ncol = K)
    Zi.path[[iter]] <- list()
    D.path[[iter]] <- list()
    Sigma.inv.path[[iter]] <- list()
    
    ## for each subject, sample beta_{i,k} from a K-dim Gaussian distribution
    beta.path[[iter]] <- foreach(i=1:N.subj, .combine='rbind', .multicombine=TRUE) %dopar% {
      beta.mean1 <- inv(t(Uk.list[[i]]) %*% t(H.list[[i]]) %*% Sigma.inv.path[[iter-1]][[i]] %*% H.list[[i]] %*% Uk.list[[i]] + 
                          Lambda.inv.path[[iter-1]]) 
      beta.mean2 <- t(Uk.list[[i]]) %*% t(H.list[[i]]) %*% Sigma.inv.path[[iter-1]][[i]] %*% data$x.demean[[i]] 
      beta.mean <- beta.mean1 %*% beta.mean2
      beta.variance <- beta.mean1
      MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
    }
    
    ## sample Lambda^-1 from a Wishart distribution
    Lambda.inv.df <- nu + N.subj + 1
    Lambda.inv.p2 <- matrix(0, nrow = K, ncol = K)
    for (i in 1:N.subj){
      Lambda.inv.p2 <- Lambda.inv.p2 + 
        (beta.path[[iter]][i,]) %*% t(beta.path[[iter]][i,])
    }
    Lambda.inv.mat <- inv(Lk + Lambda.inv.p2)  
    Lambda.inv.path[[iter]] <- MCMCpack::rwish(v = Lambda.inv.df, S = Lambda.inv.mat)
    
    ## sample Zi from a T-dim Gaussian distribution
    
    Zi.path[[iter]] <- list()
    Zi.path[[iter]] <- foreach(i=1:N.subj, .combine='c', .multicombine=TRUE) %dopar% {
      this.N.tps <- length(data$pp[[i]])
      Ai.temp <- diag(D.path[[iter-1]][[i]]) %*% Sigma.inv.path[[iter-1]][[i]]
      Ai <- diag(this.N.tps) + Ai.temp %*% diag(D.path[[iter-1]][[i]])
      ai <- Ai.temp %*% (data$x.demean[[i]]  - H.list[[i]] %*% Uk.list[[i]] %*% beta.path[[iter]][i,])
      Zi.mean <- inv(Ai) %*% ai
      Zi.varcov <- inv(Ai)
      res <- list(tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov),
                                   mean = as.numeric(Zi.mean), sigma = Zi.varcov,
                                   lower = rep(0, this.N.tps), upper = rep(Inf, this.N.tps))$samp)
      names(res) <- names(data$pp)[[i]]
      res
    }
    
    ## sample vec(D) from a T-dim Gaussian distribution
    D.path[[iter]] <- list()
    for (i in 1:N.subj){
      this.N.tps <- length(data$pp[[i]])
      B <- diag(x = 1/gamma, nrow = this.N.tps) + 
        diag(Zi.path[[iter]][[i]]) %*% Sigma.inv.path[[iter-1]][[i]] %*% diag(Zi.path[[iter]][[i]])
      b <- diag(Zi.path[[iter]][[i]]) %*% Sigma.inv.path[[iter-1]][[i]] %*%
        (data$x.demean[[i]] - H.list[[i]] %*% Uk.list[[i]] %*% beta.path[[iter]][i,]) %>% as.numeric()
      vecD.mean <- inv(B) %*% b
      vecD.varcov <- inv(B)
      D.path[[iter]][[i]] <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
    }
    
    ## sample Sigma^(-1) from a Wishart distribution
    Sigma.inv.path[[iter]] <- list()
    for (i in 1:N.subj){
      Sigma.inv.df <- twor + N.subj
      temp <- data$x.demean[[i]] - H.list[[i]] %*% Uk.list[[i]] %*% beta.path[[iter]][i,] - diag(D.path[[iter]][[i]]) %*% Zi.path[[iter]][[i]]
      Sigma.inv.mat <- inv(inv(2*kappa.list[[i]]) + temp %*% t(temp))
      Sigma.inv.path[[iter]][[i]] <- MCMCpack::rwish(v = Sigma.inv.df, S = Sigma.inv.mat)
    }
    
    #### end
    
    print(paste0("iteration ", iter, " is done"))
    
    setTxtProgressBar(pb, iter)
    
    
  }
  close(pb)
  toc()
  
  results <- list(beta = beta.path,          
                  Lambda.inv = Lambda.inv.path,            
                  Zi = Zi.path, 
                  D = D.path,
                  Sigma.inv = Sigma.inv.path,
                  H = H.list, Uk = Uk.list, Lk = Lk.list)  
  
  return(results)
  
  
}

check.traceplot.RBFPCA.sparse <- function(data, RBFPCA.sparse.result, iter.burnin, n.burnin = 1, n.iter, 
                                          i.plot, k1.plot, k2.plot){
  
  # get inputs
  beta.path <- RBFPCA.sparse.result$beta
  Lambda.inv.path <- RBFPCA.sparse.result$Lambda.inv
  Zi.path <- RBFPCA.sparse.result$Zi
  D.path <- RBFPCA.sparse.result$D
  Sigma.inv.path <- RBFPCA.sparse.result$Sigma.inv
  
  t.plot <- sample(1:length(data$pp[[i.plot]]), 1)
  
  # Get trace plots -----------------------------
  
  par(mfrow = c(2, 3))
  
  ## for beta_{i,k}
  beta.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    beta.df[iter] <- beta.path[[iter]][i.plot,k1.plot]
  }
  plot(beta.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "beta",
       main = paste0("Evaluated at i = ", i.plot, " and k = ", k1.plot))
  abline(v = iter.burnin, col = "blue")
  
  ## for Lambda^-1
  Lambda.inv.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    Lambda.inv.df[iter] <- Lambda.inv.path[[iter]][k1.plot, k2.plot]
  }
  plot(Lambda.inv.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "Lambda.inv",
       main = paste0("Evaluated at k1 = ", k1.plot, " and k2 = ", k2.plot))
  abline(v = iter.burnin, col = "blue")
  
  ## for Zi
  Zi.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    Zi.df[iter] <- Zi.path[[iter]][[i.plot]][t.plot]
  }
  plot(Zi.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "Zi",
       main = paste0("Evaluated at i = ", i.plot, " and t = ", data$pp[[i.plot]][t.plot]))
  abline(v = iter.burnin, col = "blue")
  
  ## for D
  D.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    D.df[iter] <- D.path[[iter]][[i.plot]][t.plot]
  }
  plot(D.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "D",
       main = paste0("Evaluated at t = ", data$pp[[i.plot]][t.plot]))
  abline(v = iter.burnin, col = "blue")
  
  ## for Sigma^-1
  Sigma.inv.df <- c(rep(NA, n.iter))
  for (iter in 1:n.iter){
    Sigma.inv.df[iter] <- Sigma.inv.path[[iter]][[i.plot]][t.plot, t.plot]
  }
  plot(Sigma.inv.df[-(1:n.burnin)], type = "l",
       xlab = "iteration", ylab = "Sigma.inv",
       main = paste0("Evaluated at t = ", data$pp[[i.plot]][t.plot]))
  abline(v = iter.burnin, col = "blue")
  
  par(mfrow = c(1, 1))
  
}


RBFPCA.marginal.likelihood.est <- function(RBFPCA.result, data, iter.burnin, n.iter, K.FPCA, seed){
  
  ##### First Gibbs Run #####
  
  cat("Start")
  
  beta  <- RBFPCA.res$beta[(iter.burnin+1):n.iter]
  Lambda.inv <- RBFPCA.res$Lambda.inv[(iter.burnin+1):n.iter]
  Zi  <- RBFPCA.res$Zi[(iter.burnin+1):n.iter]
  D <- RBFPCA.res$D[(iter.burnin+1):n.iter]
  Sigma.inv  <- RBFPCA.res$Sigma.inv[(iter.burnin+1):n.iter]
  
  N  <- n.iter - iter.burnin
  
  ts <- which.max(sapply(1:N, function(i) logL.fun(data, beta[[i]], Lambda.inv[[i]], Zi[[i]], D[[i]], Sigma.inv[[i]])))
  
  beta.star  <- beta[[ts]]
  Lambda.inv.star <- Lambda.inv[[ts]]
  Zi.star  <- Zi[[ts]]
  D.star <- D[[ts]]
  Sigma.inv.star <- Sigma.inv[[ts]]
  
  v1 <- logL.fun(data, beta.star, Lambda.inv.star, Zi.star, D.star, Sigma.inv.star)
  v2 <- log.prior.fun(data = data, beta = beta.star, Lambda.inv = Lambda.inv.star, 
                      Zi = Zi.star, D = D.star, Sigma.inv = Sigma.inv.star)
  chib1 <- v1 + v2
  if (chib1 == Inf | is.na(chib1)) stop(print(c(v1, v2)))
  
  
  chib2 <- chib1 - log(mean(sapply(1:(n.iter - iter.burnin), function(t) pi.beta(data, beta = beta.star, Lambda.inv = Lambda.inv[[t]],
                                                                                 Sigma.inv = Sigma.inv[[t]], Uk = RBFPCA.res$Uk, H = RBFPCA.res$H))))
  cat("First Gibbs run finished.")
  if (chib2 == Inf | is.na(chib2)) stop(print(chib1))
  
  ##### Second Gibbs Run #####
  cat("Second Gibbs run (reduced):\n")
  
  N <- n.iter
  
  BFPCA.robust.beta.fix.res <- BFPCA.robust.MC.beta.fix(data = data, num.basis = 15, num.PC = K.FPCA,
                                                        n.iter = n.iter, seed = seed, cov.index = cov.index.prior,
                                                        beta = beta.star)
  
  Lambda.inv <- BFPCA.robust.beta.fix.res$Lambda.inv[(iter.burnin+1):n.iter]
  Zi <- BFPCA.robust.beta.fix.res$Zi[(iter.burnin+1):n.iter]
  D <- BFPCA.robust.beta.fix.res$D[(iter.burnin+1):n.iter]
  Sigma.inv  <- BFPCA.robust.beta.fix.res$Sigma.inv[(iter.burnin+1):n.iter]
  
  chib3 <- chib2 - mean(sapply(1:(n.iter - iter.burnin), function(t) pi.D(data, beta = beta.star, Lambda.inv = Lambda.inv[[t]], Zi = Zi[[t]],
                                                                          D = D.star, Sigma.inv = Sigma.inv[[t]], 
                                                                          Uk = BFPCA.robust.beta.fix.res$Uk, H = BFPCA.robust.beta.fix.res$H))) 
  cat("Second Gibbs run finished.")
  if (chib3 == Inf | is.na(chib3)) stop(print(chib2))
  
  ##### Third Gibbs Run #####
  
  cat("Third Gibbs run (reduced):\n")
  
  BFPCA.robust.beta.D.fix.res <- BFPCA.robust.MC.beta.D.fix(data = data, num.basis = 15, num.PC = K.FPCA,
                                                            n.iter = n.iter, seed = seed, cov.index = cov.index.prior,
                                                            beta = beta.star, D = D.star)
  
  Lambda.inv <- BFPCA.robust.beta.D.fix.res$Lambda.inv[(iter.burnin+1):n.iter]
  Sigma.inv <- BFPCA.robust.beta.D.fix.res$Sigma.inv[(iter.burnin+1):n.iter]
  Zi <- BFPCA.robust.beta.D.fix.res$Zi[(iter.burnin+1):n.iter]
  
  chib4 <- chib3 - mean(sapply(1:(n.iter - iter.burnin), function(t) pi.Sigma.inv(data, beta = beta.star, Zi = Zi[[t]], Sigma.inv = Sigma.inv.star, D = D.star,
                                                                                  Uk = BFPCA.robust.beta.D.fix.res$Uk, H = BFPCA.robust.beta.D.fix.res$H)))
  cat("Third Gibbs run finished.")
  if (chib4 == Inf | is.na(chib4)) stop(print(chib3))
  
  ##### Fourth Gibbs Run #####
  cat("Fourth Gibbs run (reduced):\n")
  
  BFPCA.robust.beta.D.Sigma.inv.fix.res <- BFPCA.robust.MC.beta.D.Sigma.inv.fix(data = data, num.basis = 15, num.PC = K.FPCA,
                                                                                n.iter = n.iter, seed = seed, cov.index = cov.index.prior,
                                                                                beta = beta.star, D = D.star, Sigma.inv = Sigma.inv.star)
  
  Zi <- BFPCA.robust.beta.D.Sigma.inv.fix.res$Zi[(iter.burnin+1):n.iter]
  Lambda.inv <- BFPCA.robust.beta.D.Sigma.inv.fix.res$Lambda.inv[(iter.burnin+1):n.iter]
  
  chib5 <- chib4 - mean(sapply(1:(n.iter - iter.burnin), function(t) pi.Zi(data, beta = beta.star, Zi = Zi.star,
                                                                           D = D.star, Sigma.inv = Sigma.inv[[t]],
                                                                           Uk = BFPCA.robust.beta.D.Sigma.inv.fix.res$Uk, 
                                                                           H = BFPCA.robust.beta.D.Sigma.inv.fix.res$H)))
  cat("Fourth Gibbs run finished.")
  if (chib5 == Inf | is.na(chib5)) stop(print(chib4))
  
  return(chib5)
  
}

RBFPCA.Bayes.factor <- function(RBFPCA.result1, K.FPCA1, RBFPCA.result2, K.FPCA2, data, iter.burnin, n.iter, seed){
  val1 <- RBFPCA.marginal.likelihood.est(RBFPCA.res1, data, iter.burnin, n.iter, K.FPCA1, seed)
  val2 <- RBFPCA.marginal.likelihood.est(RBFPCA.res2, data, iter.burnin, n.iter, K.FPCA2, seed)
  return(exp(val1-val2))
}