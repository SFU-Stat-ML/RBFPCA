
#####################################################################

RBFPCA.ASMC.mvn <- function(data, num.basis, num.PC, tuning_param, seed, cov.index, 
                            N.particles, n.cores) {
  # data: Input data matrix of dimension N.tps * N.subj
  # num.basis: Number of basis functions (J)
  # num.PC: Number of principal components (K)
  # tuning_param: List of tuning parameters for the ASMC algorithm: 1) eps - resampling threshold,
  #  2) phi - target rCESS for adaptive resampling scheme, 3) alhpa - annealing scheme. Only supply one of phi or alpha, not both.
  # seed: Seed for random number generation
  # cov.index: Covariance function index
  # N.particles: Number of particles for the SMC algorithm
  # n.cores: Number of cores for the SMC algorithm
  
  set.seed(seed)
  
  # Get inputs -----------------------------
  N.subj <- ncol(data) # Number of subjects
  N.tps <- nrow(data)  # Number of time points
  
  J <- num.basis  # Number of basis functions
  K <- num.PC     # Number of principal components
  
  # Define the hyperparameters -----------------------------
  CovFun <- function(s, t, cov.index) {
    if (cov.index == 1) {
      exp(-3 * (t - s)^2)    # Exponential covariance function
    } else {
      min(s + 1, t + 1)     # Linear covariance function
    }
  }
  
  ## U_k, L_k ##
  # Construct covariance matrix Omega.star
  timepts <- seq(-1, 1, length.out = N.tps)
  Omega.star <- matrix(nrow = N.tps, ncol = N.tps)
  if (cov.index == 0) {
    Omega.star <- cov(t(data))   # Empirical covariance matrix
  } else {
    for (i in 1:N.tps) {
      for (j in 1:N.tps) {
        Omega.star[i, j] <- CovFun(s = timepts[i], t = timepts[j], cov.index)
      }
    }
  }
  
  # Compute the basis matrix H using Legendre polynomials
  normalized.p.list <- legendre.polynomials(J - 1, normalized = TRUE)
  H <- matrix(nrow = N.tps, ncol = J)
  for (i in 1:length(normalized.p.list)) {
    H[, i] <- predict(normalized.p.list[[i]], newdata = timepts)
  }
  HH.inv <- inv(t(H) %*% H)
  
  # Compute eigen decomposition for Uk and Lk
  xi <- HH.inv %*% t(H) %*% Omega.star %*% H %*% HH.inv
  Uk <- eigen(xi)$vectors[, 1:K]
  Lk <- diag(eigen(xi)$values[1:K])
  HUk = H %*% Uk
  UHHU <- t(HUk) %*% HUk
  Lk.inv <- solve(Lk)
  
  # Set hyperparameters for the Gibbs sampler
  
  ## nu ##
  nu <- 2 * K  # K is Number of principal components  
  
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
  kappa <- (100 / twor) * inv(R)
  
  ## gamma ##
  gamma <- 10
  
  prior <- function(theta) {
    
    # Precompute values outside the loop
    log_prior <- log(MCMCpack::dwish(theta$Omega.inv, v = nu, S = Lk.inv))
    log_prior <- log_prior + LaplacesDemon::dwishart(Omega = theta$Sigma.inv, nu = twor, S = 2 * kappa, log = TRUE)
    
    Omega <- solve(theta$Omega.inv)  # Precompute the inverse once
    mean_zero <- rep(0, K)  # Precompute repeated vectors
    mu_zi <- rep(0, N.tps)
    diag_sigma <- diag(1, N.tps)
    lb <- rep(0, N.tps)
    ub <- rep(Inf, N.tps)
    
    # Use vectorized operations where possible
    log_prior <- log_prior + sum(mvtnorm::dmvnorm(theta$beta, mean_zero, Omega, log = TRUE))
    log_prior <- log_prior + sum(TruncatedNormal::dtmvnorm(theta$Zi, mu = mu_zi, sigma = diag_sigma, lb = lb, ub = ub, log = TRUE))
    
    # Sum dnorm for theta$D outside the loop
    log_prior <- log_prior + sum(dnorm(theta$D, mean = 0, sd = 10, log = TRUE))
    
    return(log_prior)
  }
  
  likelihood <- function(theta){
    log_likelihood <- 0
    sigma_inv <- theta$Sigma.inv
    Sigma <- solve(sigma_inv)
    d <- diag(theta$D)
    for (i in 1:N.subj) {
      mean_part <- HUk %*% theta$beta[i, ]
      zi <- theta$Zi[i, ]
      log_likelihood <- log_likelihood + mvtnorm::dmvnorm(data[, i], mean = mean_part + d %*% zi, sigma = Sigma, log = TRUE)
    }
    return(log_likelihood)
  }
  
  
  # Define initial particles -----------------------------
  
  # initialize storage list
  alpha <- list() # a list for tempering parameters
  ESS <- list()   # a list for ESS
  logZ <- list()  # list for log-normalizing constant
  W <- list()
  
  # check if adaptive annealing scheme is being used
  adaptive_alpha <- "phi" %in% names(tuning_param)
  
  # initialize values
  r <- 1                  # SMC iteration
  logZ[[r]] <- 0
  
  if(adaptive_alpha){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(tuning_param$alpha)
  }
  
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl,{
    library(MCMCpack)
    library(tmvnsim)
    library(mvtnorm)
    library(mnormt)
    library(polynom)
    library(MASS)
    library(TruncatedNormal)
    library(LaplacesDemon)
  })
  
  # initialize particles
  init <- function(j)
  {
    set.seed(j)
    init_particle <- list(
      beta = matrix(0, nrow = N.subj, ncol = K),
      Omega.inv = MCMCpack::rwish(v = nu, S = Lk.inv),
      Zi = matrix(data = NA, nrow = N.subj, ncol = N.tps),
      D = rnorm(n = N.tps, mean = 0, sd = 10),
      Sigma.inv = MCMCpack::rwish(v = twor, S = 2*kappa)
    )
    
    Omega <- solve(init_particle$Omega.inv)
    for (i in 1:N.subj) {
      init_particle$beta[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega)
      init_particle$Zi[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps, mean = rep(0,N.tps), sigma = diag(1,N.tps), lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
    }
    return(init_particle)
  }
  previous_particles <- parLapply(cl, 1:N.particles, function(k) init(k))
  
  ## Initialize weights
  # normalized weights
  # Particles were sampled from the same reference distribution. They have the equal weights.
  W[[r]] <- rep(1/N.particles,N.particles)
  logW <- log(W[[r]])
  
  # unnormalized weights
  w <- rep(1,N.particles)
  logw <- rep(0,N.particles)
  
  # Annealed SMC Algorithm -----------------------------
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    # evaluate the log-likelihood for updating alpha
    u <- rep(0, N.particles)   # incremental log-importance weights
    u <- sapply(1:N.particles, function(k){
      logL <- likelihood(previous_particles[[k]])
      return(logL)
    })
    
    if(adaptive_alpha){
      # update alpha with bisection
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$phi )
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n")
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]]
    }
    
    Gibbs_HM_Move <- function(theta)
    {
      theta_new <- theta
      # Randomly pick one parameter to update
      chosen_param <- sample(names(theta_new), 1)
      
      proposal_function_beta <- function() {
        new_beta <- theta$beta
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        beta.variance <- solve(t(HUk) %*% theta$Sigma.inv %*% HUk + theta$Omega.inv)
        beta.mean <- beta.variance %*% t(HUk) %*% theta$Sigma.inv %*% data[, i]
        new_beta[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
        return(new_beta)
      }
      
      proposal_function_Omega_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        Omega.inv.df <- nu + N.subj + 1
        Omega.inv.p2 <- matrix(0, nrow = K, ncol = K)
        for (i in 1:N.subj) {
          Omega.inv.p2 <- Omega.inv.p2 + (theta$beta[i, ]) %*% t(theta$beta[i, ])
        }
        Omega.inv.mat <- solve(Lk + Omega.inv.p2)
        new.Omega.inv <- MCMCpack::rwish(v = Omega.inv.df, S = Omega.inv.mat)
        return(new.Omega.inv)
      }
      
      proposal_function_Zi <- function() {
        new.Zi <- theta$Zi
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        Ai.temp <- diag(theta$D) %*% theta$Sigma.inv
        Ai <- diag(N.tps) + Ai.temp %*% diag(theta$D)
        ai <- Ai.temp %*% (data[, i] - HUk %*% theta$beta[i, ])
        Zi.varcov <- solve(Ai)
        Zi.mean <- Zi.varcov %*% ai
        new.Zi[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), mean = Zi.mean, sigma = Zi.varcov, lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
        return(new.Zi)
      }
      
      proposal_function_D <- function() {
        
        B.temp <- matrix(data = 0, nrow = N.tps, ncol = N.tps) + diag(x = 1/gamma, nrow = N.tps)
        b.temp <- matrix(data = 0, nrow = N.tps, ncol = 1)
        for (i in 1:N.subj){
          B.temp <- B.temp + diag(theta$Zi[i, ]) %*% theta$Sigma.inv %*% diag(theta$Zi[i, ])
          b.temp <- b.temp + diag(theta$Zi[i, ]) %*% theta$Sigma.inv %*%
            (data[, i] - HUk %*% theta$beta[i,])
        }
        vecD.mean <- solve(B.temp) %*% b.temp
        vecD.varcov <- solve(B.temp)
        new.D <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
        return(new.D)
      }
      
      proposal_function_Sigma_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        Sigma.inv.df <- twor + N.subj
        Sigma.inv.mat <- solve(2 * kappa)
        for (i in 1:N.subj) {
          temp <- data[, i] - HUk %*% theta$beta[i, ] - diag(theta$D) %*% theta$Zi[i, ]
          Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
        }
        new.Sigma.inv <- MCMCpack::rwish(v = Sigma.inv.df, S = solve(Sigma.inv.mat))
        return(new.Sigma.inv)
      }
      
      # Map each parameter to its corresponding proposal function
      proposal_functions <- list(
        beta = proposal_function_beta,
        Omega.inv = proposal_function_Omega_inv,
        Zi = proposal_function_Zi,
        D =  proposal_function_D,
        Sigma.inv = proposal_function_Sigma_inv
      )
      
      # Update the chosen parameter using its respective proposal function
      theta_new[[chosen_param]] <- proposal_functions[[chosen_param]]()
      
      # compute numerator and denominator of the MH ratio
      num <- (likelihood(theta_new)) * alpha[[r]] +  prior(theta_new)
      den <- (likelihood(theta)) * alpha[[r]] + prior(theta)
      
      # compute ratio
      ratio <- min(1, exp(num - den))
      
      # accept/reject step
      if(runif(1) < ratio){
        theta = theta_new
      } else{
        theta = theta
      }
      
      return(theta)
    }
  
    particles <- parLapply(cl, previous_particles, Gibbs_HM_Move)
    
    # compute the ESS
    log_incremental_w <- alphaDiff * u
    logw <- log_incremental_w + logw  # log unnormalized weights
    logmax <- max(logw)
    logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  
    W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
    logW <- log(W[[r]])                                # log normalized weights
    cat("logZ:",logZ[[r]],"\n")
    ESS[[r]] <- rESS(logW)
    
    # resample if ESS below threshold
    if( ESS[[r]]<tuning_param$eps )
    {
      cat("Resample: ESS=", ESS[[r]], '\n')
      ancestors <- systematic_resample( W[[r]] )
      particles <-  particles[ancestors]
      W[[r]] <- rep(1/N.particles,N.particles)
      logW <- log(W[[r]])
      w <- rep(1,N.particles)
      logw <- rep(0,N.particles)
    }
    previous_particles <- particles
  }
 
  stopCluster(cl)
  
  return(list(particles = particles,
              alpha = alpha,
              ESS = ESS,
              logZ = logZ,
              W = W,
              H = H, Uk = Uk, Lk = Lk))
}


RBFPCA.ASMC.mvt <- function(data, num.basis, num.PC, tuning_param, seed, cov.index, 
                            N.particles, n.cores) {
  # data: Input data matrix of dimension N.tps * N.subj
  # num.basis: Number of basis functions (J)
  # num.PC: Number of principal components (K)
  # tuning_param: List of tuning parameters for the ASMC algorithm: 1) eps - resampling threshold,
  #  2) phi - target rCESS for adaptive resampling scheme, 3) alhpa - annealing scheme. Only supply one of phi or alpha, not both.
  # seed: Seed for random number generation
  # cov.index: Covariance function index
  # N.particles: Number of particles for the SMC algorithm
  # n.cores: Number of cores for the SMC algorithm
  
  set.seed(seed)
  
  # Get inputs -----------------------------
  N.subj <- ncol(data) # Number of subjects
  N.tps <- nrow(data)  # Number of time points
  
  J <- num.basis  # Number of basis functions
  K <- num.PC     # Number of principal components
  
  # Define the hyperparameters -----------------------------
  CovFun <- function(s, t, cov.index) {
    if (cov.index == 1) {
      exp(-3 * (t - s)^2)    # Exponential covariance function
    } else {
      min(s + 1, t + 1)     # Linear covariance function
    }
  }
  
  ## U_k, L_k ##
  # Construct covariance matrix Omega.star
  timepts <- seq(-1, 1, length.out = N.tps)
  Omega.star <- matrix(nrow = N.tps, ncol = N.tps)
  if (cov.index == 0) {
    Omega.star <- cov(t(data))   # Empirical covariance matrix
  } else {
    for (i in 1:N.tps) {
      for (j in 1:N.tps) {
        Omega.star[i, j] <- CovFun(s = timepts[i], t = timepts[j], cov.index)
      }
    }
  }
  
  # Compute the basis matrix H using Legendre polynomials
  normalized.p.list <- legendre.polynomials(J - 1, normalized = TRUE)
  H <- matrix(nrow = N.tps, ncol = J)
  for (i in 1:length(normalized.p.list)) {
    H[, i] <- predict(normalized.p.list[[i]], newdata = timepts)
  }
  HH.inv <- inv(t(H) %*% H)
  
  # Compute eigen decomposition for Uk and Lk
  xi <- HH.inv %*% t(H) %*% Omega.star %*% H %*% HH.inv
  Uk <- eigen(xi)$vectors[, 1:K]
  Lk <- diag(eigen(xi)$values[1:K])
  HUk = H %*% Uk
  UHHU <- t(HUk) %*% HUk
  Lk.inv <- solve(Lk)
  
  # Set hyperparameters for the Gibbs sampler
  
  ## nu ##
  nu <- 2 * K  # K is Number of principal components  
  
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
  kappa <- (100 / twor) * inv(R)
  
  ## gamma ##
  gamma <- 10
  
  prior <- function(theta) {
    
    # Precompute values outside the loop
    log_prior <- log(MCMCpack::dwish(theta$Omega.inv, v = nu, S = Lk.inv))
    log_prior <- log_prior + LaplacesDemon::dwishart(Omega = theta$Sigma.inv, nu = twor, S = 2 * kappa, log = TRUE)
    
    Omega <- solve(theta$Omega.inv)  # Precompute the inverse once
    mean_zero <- rep(0, K)  # Precompute repeated vectors
    mu_zi <- rep(0, N.tps)
    diag_sigma <- diag(1, N.tps)
    lb <- rep(0, N.tps)
    ub <- rep(Inf, N.tps)
    
    # Use vectorized operations where possible
    log_prior <- log_prior + sum(mvtnorm::dmvnorm(theta$beta, mean_zero, Omega, log = TRUE))
    log_prior <- log_prior + sum(TruncatedNormal::dtmvnorm(theta$Zi, mu = mu_zi, sigma = diag_sigma, lb = lb, ub = ub, log = TRUE))
    
    # Sum dnorm for theta$D outside the loop
    log_prior <- log_prior + sum(dnorm(theta$D, mean = 0, sd = gamma, log = TRUE))
    
    # Sum v and wi
    log_prior <- log_prior + sum(LaplacesDemon::dtrunc(theta$v, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1))
    log_prior <- log_prior + sum(dgamma(theta$wi, shape = theta$v/2, rate = theta$v/2))
    
    return(log_prior)
  }
  
  likelihood <- function(theta){
    log_likelihood <- 0
    sigma_inv <- theta$Sigma.inv
    Sigma <- solve(sigma_inv)
    d <- diag(theta$D)
    wi <- theta$wi
    for (i in 1:N.subj) {
      mean_part <- HUk %*% theta$beta[i, ]
      zi <- theta$Zi[i, ]
      log_likelihood <- log_likelihood + mvtnorm::dmvnorm(data[, i], mean = mean_part + d %*% zi, sigma = Sigma/wi[i], log = TRUE)
    }
    return(log_likelihood)
  }
  
  
  # Define initial particles -----------------------------
  
  # initialize storage list
  alpha <- list() # a list for tempering parameters
  ESS <- list()   # a list for ESS
  logZ <- list()  # list for log-normalizing constant
  W <- list()
  
  # check if adaptive annealing scheme is being used
  adaptive_alpha <- "phi" %in% names(tuning_param)
  
  # initialize values
  r <- 1                  # SMC iteration
  logZ[[r]] <- 0
  
  if(adaptive_alpha){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(tuning_param$alpha)
  }
  
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl,{
    library(MCMCpack)
    library(tmvnsim)
    library(mvtnorm)
    library(mnormt)
    library(polynom)
    library(MASS)
    library(TruncatedNormal)
    library(LaplacesDemon)
  })
  
  # initialize particles
  init <- function(j)
  {
    set.seed(j)
    wi_func <- function(x){
      rgamma(1,shape = x/2, rate = x/2)
    }
    init_particle <- list(
      beta = matrix(0, nrow = N.subj, ncol = K),
      Omega.inv = MCMCpack::rwish(v = nu, S = Lk.inv), 
      Zi = matrix(data = NA, nrow = N.subj, ncol = N.tps),
      D = rnorm(n = N.tps, mean = 0, sd = 10),
      Sigma.inv = MCMCpack::rwish(v = twor, S = 2*kappa)
    )
    
    init_particle$vi = LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1)
    init_particle$wi <- sapply(init_particle$vi, wi_func)
    
    Omega <- solve(init_particle$Omega.inv)
    for (i in 1:N.subj) {
      init_particle$beta[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega)
      init_particle$Zi[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps, mean = rep(0,N.tps), sigma = diag(1,N.tps), lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
    }
    return(init_particle)
  }
  
  previous_particles <- parLapply(cl, 1:N.particles, function(k) init(k))
  
  ## Initialize weights
  # normalized weights
  # Particles were sampled from the same reference distribution. They have the equal weights.
  W[[r]] <- rep(1/N.particles,N.particles)
  logW <- log(W[[r]])
  
  # unnormalized weights
  w <- rep(1,N.particles)
  logw <- rep(0,N.particles)
  
  # Annealed SMC Algorithm -----------------------------
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    # evaluate the log-likelihood for updating alpha
    u <- rep(0, N.particles)   # incremental log-importance weights
    u <- sapply(1:N.particles, function(k){
      logL <- likelihood(previous_particles[[k]])
      return(logL)
    })
    
    if(adaptive_alpha){
      # update alpha with bisection
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$phi )
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n")
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]]
    }
    
    Gibbs_HM_Move <- function(theta)
    {

      theta_new <- theta
      # Randomly pick one parameter to update
      param_names <- names(theta_new)[!names(theta_new) == "vi"]
      chosen_param <- sample(param_names, 1)
      
      proposal_function_beta <- function() {
        new_beta <- theta$beta
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        beta.variance <- solve(t(HUk) %*% theta$Sigma.inv %*% HUk + theta$Omega.inv)
        beta.mean <- beta.variance %*% t(HUk) %*% theta$Sigma.inv %*% data[, i]
        new_beta[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
        return(new_beta)
      }
      
      proposal_function_Omega_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        Omega.inv.df <- nu + N.subj + 1
        Omega.inv.p2 <- matrix(0, nrow = K, ncol = K)
        for (i in 1:N.subj) {
          Omega.inv.p2 <- Omega.inv.p2 + (theta$beta[i, ]) %*% t(theta$beta[i, ])
        }
        Omega.inv.mat <- solve(Lk + Omega.inv.p2)
        new.Omega.inv <- MCMCpack::rwish(v = Omega.inv.df, S = Omega.inv.mat)
        return(new.Omega.inv)
      }
      
      proposal_function_Zi <- function() {
        new.Zi <- theta$Zi
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        Ai.temp <- diag(theta$D) %*% theta$Sigma.inv
        Ai <- diag(N.tps) + Ai.temp %*% diag(theta$D)
        ai <- Ai.temp %*% (data[, i] - HUk %*% theta$beta[i, ])
        Zi.varcov <- solve(Ai)
        Zi.mean <- Zi.varcov %*% ai
        new.Zi[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), mean = Zi.mean, sigma = Zi.varcov, lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
        return(new.Zi)
      }
      
      proposal_function_D <- function() {
        add <- function(x) Reduce("+", x)
        BB <- matrix(data = NA, nrow = N.tps * N.subj, ncol = N.tps)
        bb <- matrix(data = NA, nrow = N.subj, ncol = N.tps)
        for (i in 1:N.subj) {
          B.temp <- diag(theta$Zi[i, ]) %*% theta$Sigma.inv %*% diag(theta$Zi[i, ])
          b.temp <- diag(theta$Zi[i, ]) %*% theta$Sigma.inv %*% (data[, i] - HUk %*% theta$beta[i, ])
          BB[((i - 1) * N.tps + 1):(i * N.tps), ] <- B.temp
          bb[i, ] <- b.temp
        }
        res <- list(B.temp = BB, b.temp = bb)
        
        number_of_chunks <- nrow(res$B.temp) / N.tps
        B.temp.list <- lapply(split(res$B.temp, rep(1:number_of_chunks, each = NROW(res$B.temp) / number_of_chunks)), function(a) matrix(a, ncol = NCOL(res$B.temp)))
        B.temp <- add(B.temp.list) + diag(x = 1 / gamma, nrow = N.tps)
        b.temp <- colSums(res$b.temp)
        vecD.varcov <- solve(B.temp)
        vecD.mean <- vecD.varcov %*% b.temp
        
        new.D <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
        return(new.D)
      }
      
      proposal_function_Sigma_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        Sigma.inv.df <- twor + N.subj
        Sigma.inv.mat <- solve(2 * kappa)
        for (i in 1:N.subj) {
          temp <- data[, i] - HUk %*% theta$beta[i, ] - diag(theta$D) %*% theta$Zi[i, ]
          Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
        }
        new.Sigma.inv <- MCMCpack::rwish(v = Sigma.inv.df, S = solve(Sigma.inv.mat))
        return(new.Sigma.inv)
      }
      
      proposal_function_wi <- function() {
        new.wi <- theta$wi
        ## sample wi from a Gamma distribution
        vi <- LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1)
        wi.temp <- apply(data, 2, var)
        for (i in 1:N.subj){
          wi.shape <- vi[i]/2 + N.subj/2
          wi.rate <- vi[i]/2 +  N.subj * wi.temp[i]/ 2
          new.wi[i] <- rgamma(1, shape = wi.shape, rate = wi.rate)
        }
        return(new.wi)
      }
      
      # Map each parameter to its corresponding proposal function
      proposal_functions <- list(
        beta = proposal_function_beta,
        Omega.inv = proposal_function_Omega_inv,
        Zi = proposal_function_Zi,
        D =  proposal_function_D,
        Sigma.inv = proposal_function_Sigma_inv,
        wi = proposal_function_wi
      )
      
      # Update the chosen parameter using its respective proposal function
      theta_new[[chosen_param]] <- proposal_functions[[chosen_param]]()
      
      # compute numerator and denominator of the MH ratio
      num <- (likelihood(theta_new)) * alpha[[r]] +  prior(theta_new)
      den <- (likelihood(theta)) * alpha[[r]] + prior(theta)
      
      # compute ratio
      ratio <- min(1, exp(num - den))
      
      # accept/reject step
      if(runif(1) < ratio){
        theta = theta_new
      } else{
        theta = theta
      }
      
      return(theta)
    }
    
    particles <- parLapply(cl, previous_particles, Gibbs_HM_Move)
    
    # compute the ESS
    log_incremental_w <- alphaDiff * u
    logw <- log_incremental_w + logw  # log unnormalized weights
    logmax <- max(logw)
    logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  
    W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
    logW <- log(W[[r]])                                # log normalized weights
    cat("logZ:",logZ[[r]],"\n")
    ESS[[r]] <- rESS(logW)
    
    # resample if ESS below threshold
    if( ESS[[r]]<tuning_param$eps )
    {
      cat("Resample: ESS=", ESS[[r]], '\n')
      ancestors <- systematic_resample( W[[r]] )
      particles <-  particles[ancestors]
      W[[r]] <- rep(1/N.particles,N.particles)
      logW <- log(W[[r]])
      w <- rep(1,N.particles)
      logw <- rep(0,N.particles)
    }
    previous_particles <- particles
    
  }

  stopCluster(cl)
  
  return(list(particles = particles,
              alpha = alpha,
              ESS = ESS,
              logZ = logZ,
              W = W,
              H = H, Uk = Uk, Lk = Lk))
}


RBFPCA.ASMC.mm <- function(data, num.basis, num.PC, tuning_param, seed, cov.index, 
                           N.particles, n.cores) {
  # data: Input data matrix of dimension N.tps * N.subj
  # num.basis: Number of basis functions (J)
  # num.PC: Number of principal components (K)
  # tuning_param: List of tuning parameters for the ASMC algorithm: 1) eps - resampling threshold,
  #  2) phi - target rCESS for adaptive resampling scheme, 3) alhpa - annealing scheme. Only supply one of phi or alpha, not both.
  # seed: Seed for random number generation
  # cov.index: Covariance function index
  # N.particles: Number of particles for the SMC algorithm
  # n.cores: Number of cores for the SMC algorithm
  
  set.seed(seed)
  
  # Get inputs -----------------------------
  N.subj <- ncol(data) # Number of subjects
  N.tps <- nrow(data)  # Number of time points
  
  J <- num.basis  # Number of basis functions
  K <- num.PC     # Number of principal components
  
  # Define the hyperparameters -----------------------------
  CovFun <- function(s, t, cov.index) {
    if (cov.index == 1) {
      exp(-3 * (t - s)^2)    # Exponential covariance function
    } else {
      min(s + 1, t + 1)     # Linear covariance function
    }
  }
  
  ## U_k, L_k ##
  # Construct covariance matrix Omega.star
  timepts <- seq(-1, 1, length.out = N.tps)
  Omega.star <- matrix(nrow = N.tps, ncol = N.tps)
  if (cov.index == 0) {
    Omega.star <- cov(t(data))   # Empirical covariance matrix
  } else {
    for (i in 1:N.tps) {
      for (j in 1:N.tps) {
        Omega.star[i, j] <- CovFun(s = timepts[i], t = timepts[j], cov.index)
      }
    }
  }
  
  # Compute the basis matrix H using Legendre polynomials
  normalized.p.list <- legendre.polynomials(J - 1, normalized = TRUE)
  H <- matrix(nrow = N.tps, ncol = J)
  for (i in 1:length(normalized.p.list)) {
    H[, i] <- predict(normalized.p.list[[i]], newdata = timepts)
  }
  HH.inv <- inv(t(H) %*% H)
  
  # Compute eigen decomposition for Uk and Lk
  xi <- HH.inv %*% t(H) %*% Omega.star %*% H %*% HH.inv
  Uk <- eigen(xi)$vectors[, 1:K]
  Lk <- diag(eigen(xi)$values[1:K])
  HUk = H %*% Uk
  UHHU <- t(HUk) %*% HUk
  Lk.inv <- solve(Lk)
  
  # Set hyperparameters for the Gibbs sampler
  
  ## nu ##
  nu <- 2 * K  # K is Number of principal components  
  
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
  kappa <- (100 / twor) * inv(R)
  
  ## gamma ##
  gamma <- 10
  
  prior <- function(theta) {
    
    # Precompute values outside the loop
    log_prior <- log(MCMCpack::dwish(theta$Omega1.inv, v = nu, S = Lk.inv)) +
      log(MCMCpack::dwish(theta$Omega2.inv, v = nu, S = Lk.inv))
    log_prior <- log_prior + LaplacesDemon::dwishart(Omega = theta$Sigma1.inv, nu = twor, S = 2 * kappa, log = TRUE) +
      LaplacesDemon::dwishart(Omega = theta$Sigma2.inv, nu = twor, S = 2 * kappa, log = TRUE)
    
    Omega1 <- solve(theta$Omega1.inv)  # Precompute the inverse once
    Omega2 <- solve(theta$Omega2.inv)  # Precompute the inverse once
    mean_zero <- rep(0, K)  # Precompute repeated vectors
    mu_zi <- rep(0, N.tps)
    diag_sigma <- diag(1, N.tps)
    lb <- rep(0, N.tps)
    ub <- rep(Inf, N.tps)
    
    # Use vectorized operations where possible
    log_prior <- log_prior + sum(mvtnorm::dmvnorm(theta$beta1, mean_zero, Omega1, log = TRUE)) +
      sum(mvtnorm::dmvnorm(theta$beta2, mean_zero, Omega2, log = TRUE))
    log_prior <- log_prior + sum(TruncatedNormal::dtmvnorm(theta$Zi1, mu = mu_zi, sigma = diag_sigma, lb = lb, ub = ub, log = TRUE)) +
      sum(TruncatedNormal::dtmvnorm(theta$Zi2, mu = mu_zi, sigma = diag_sigma, lb = lb, ub = ub, log = TRUE))
    
    # Sum dnorm for theta$D outside the loop
    log_prior <- log_prior + sum(dnorm(theta$D1, mean = 0, sd = 10, log = TRUE)) +
      sum(dnorm(theta$D2, mean = 0, sd = 10, log = TRUE))
    
    # Sum v and wi
    log_prior <- log_prior + sum(LaplacesDemon::dtrunc(theta$v, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1))
    log_prior <- log_prior + sum(dgamma(theta$wi, shape = theta$v/2, rate = theta$v/2))
    
    # Sum pi
    log_prior <- log_prior + sum(dbeta(theta$pi1, 1, 1))
    
    return(log_prior)
  }
  
  likelihood <- function(theta){
    log_likelihood <- 0
    sigma1_inv <- theta$Sigma1.inv
    Sigma1 <- solve(sigma1_inv)
    d1 <- diag(theta$D1)
    sigma2_inv <- theta$Sigma2.inv
    Sigma2 <- solve(sigma2_inv)
    d2 <- diag(theta$D2)
    wi <- theta$wi
    tau <- theta$tau
    pi1 <- theta$pi1
    for (i in 1:N.subj) {
      if (tau[i] == 1){
        mean_part2 <- HUk %*% theta$beta2[i, ]
        zi2 <- theta$Zi2[i, ]
        log_likelihood <- log_likelihood + log(1-pi1) + mvtnorm::dmvnorm(data[, i], mean = mean_part2 + d2 %*% zi2, sigma = Sigma2/wi[i], log = TRUE)
      } else{
        mean_part1 <- HUk %*% theta$beta1[i, ]
        zi1 <- theta$Zi1[i, ]
        log_likelihood <- log_likelihood + log(pi1) + mvtnorm::dmvnorm(data[, i], mean = mean_part1 + d1 %*% zi1, sigma = Sigma1, log = TRUE)
      }
      
    }
    return(log_likelihood)
  }
  
  
  # Define initial particles -----------------------------
  
  # initialize storage list
  alpha <- list() # a list for tempering parameters
  ESS <- list()   # a list for ESS
  logZ <- list()  # list for log-normalizing constant
  W <- list()
  
  # check if adaptive annealing scheme is being used
  adaptive_alpha <- "phi" %in% names(tuning_param)
  
  # initialize values
  r <- 1                  # SMC iteration
  logZ[[r]] <- 0
  
  if(adaptive_alpha){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(tuning_param$alpha)
  }
  
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl,{
    library(MCMCpack)
    library(tmvnsim)
    library(mvtnorm)
    library(mnormt)
    library(polynom)
    library(MASS)
    library(TruncatedNormal)
    library(LaplacesDemon)
  })
  
  # initialize particles
  init <- function(j)
  {
    set.seed(j)
    ReMap.mm <- function(x){
      x.unique <- unique(x)
      x.remap <- vector(mode = "integer", length = length(x))
      for (i in 1:length(x.unique)){
        x.remap[x == x.unique[i]] <- i
      }
      return(x.remap)
    }
    
    wi_func <- function(x){
      rgamma(1,shape = x/2, rate = x/2)
    }
    init_particle <- list(
      beta1 = matrix(0, nrow = N.subj, ncol = K),
      beta2 = matrix(0, nrow = N.subj, ncol = K),
      Omega1.inv = MCMCpack::rwish(v = nu, S = Lk.inv), 
      Omega2.inv = MCMCpack::rwish(v = nu, S = Lk.inv),
      Zi1 = matrix(data = NA, nrow = N.subj, ncol = N.tps),
      Zi2 = matrix(data = NA, nrow = N.subj, ncol = N.tps),
      D1 = rnorm(n = N.tps, mean = 0, sd = 10),
      D2 = rnorm(n = N.tps, mean = 0, sd = 10),
      Sigma1.inv = MCMCpack::rwish(v = twor, S = 2*kappa),
      Sigma2.inv = MCMCpack::rwish(v = twor, S = 2*kappa)
    )
    
    init_particle$vi = LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1)
    init_particle$wi <- sapply(init_particle$vi, wi_func)
    
    init_particle$pi1 <- rbeta(1,1,1)
    init_particle$tau <- ReMap.mm((LaplacesDemon::rbern(N.subj, init_particle$pi1) + 1))
    while (length(unique(init_particle$tau)) == 1){
      init_particle$pi1 <- rbeta(1,1,1)
      init_particle$tau <- ReMap.mm((LaplacesDemon::rbern(N.subj, init_particle$pi1) + 1))
    }
    
    Omega1 <- solve(init_particle$Omega1.inv)
    Omega2 <- solve(init_particle$Omega2.inv)
    for (i in 1:N.subj) {
      init_particle$beta1[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega1)
      init_particle$beta2[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega2)
      init_particle$Zi1[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps, mean = rep(0,N.tps), 
                                                 sigma = diag(1,N.tps), lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
      init_particle$Zi2[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps, mean = rep(0,N.tps), 
                                                 sigma = diag(1,N.tps), lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
    }
    return(init_particle)
  }
  
  previous_particles <- parLapply(cl, 1:N.particles, function(k) init(k))
  
  ## Initialize weights
  # normalized weights
  # Particles were sampled from the same reference distribution. They have the equal weights.
  W[[r]] <- rep(1/N.particles,N.particles)
  logW <- log(W[[r]])
  
  # unnormalized weights
  w <- rep(1,N.particles)
  logw <- rep(0,N.particles)
  
  # Annealed SMC Algorithm -----------------------------
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    # evaluate the log-likelihood for updating alpha
    u <- rep(0, N.particles)   # incremental log-importance weights
    u <- sapply(1:N.particles, function(k){
      logL <- likelihood(previous_particles[[k]])
      return(logL)
    })
    
    if(adaptive_alpha){
      # update alpha with bisection
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$phi )
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n")
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]]
    }
    
    Gibbs_HM_Move <- function(theta)
    {
      
      ReMap.mm <- function(x){
        x.unique <- unique(x)
        x.remap <- vector(mode = "integer", length = length(x))
        for (i in 1:length(x.unique)){
          x.remap[x == x.unique[i]] <- i
        }
        return(x.remap)
      }
      
      theta_new <- theta
      # Randomly pick one parameter to update
      param_names <- names(theta_new)[!names(theta_new) == "vi"]
      chosen_param <- sample(param_names, 1)
      
      proposal_function_beta1 <- function(update.all = FALSE) {
        new_beta1 <- theta$beta1
        tau <- theta$tau
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new_beta1 <- theta$beta1
        } else{
          if (update.all){
            for (i in subj1){
              beta1.variance <- solve(t(HUk) %*% theta$Sigma1.inv %*% HUk + theta$Omega1.inv)
              beta1.mean <- beta1.variance %*% t(HUk) %*% theta$Sigma1.inv %*% data[, i]
              new_beta1[i, ] <- MASS::mvrnorm(n = 1, mu = beta1.mean, Sigma = beta1.variance)
            }
          } else{
            # Sample one value from component 1 (mvn)
            i <- sample(subj1, 1)
            beta1.variance <- solve(t(HUk) %*% theta$Sigma1.inv %*% HUk + theta$Omega1.inv)
            beta1.mean <- beta1.variance %*% t(HUk) %*% theta$Sigma1.inv %*% data[, i]
            new_beta1[i, ] <- MASS::mvrnorm(n = 1, mu = beta1.mean, Sigma = beta1.variance)
          }
        }
        return(new_beta1)
      }
      
      proposal_function_beta2 <- function(update.all = FALSE) {
        new_beta2 <- theta$beta2
        tau <- theta$tau
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new_beta2 <- theta$beta2
        } else{
          if (update.all){
            for (i in subj2){
              beta2.variance <- solve(t(HUk) %*% theta$Sigma2.inv %*% HUk + theta$Omega2.inv)
              beta2.mean <- beta2.variance %*% t(HUk) %*% theta$Sigma2.inv %*% data[, i]
              new_beta2[i, ] <- MASS::mvrnorm(n = 1, mu = beta2.mean, Sigma = beta2.variance)
            }
          } else{
            # Sample one value from component 2 (mvt)
            i <- sample(subj2, 1)
            beta2.variance <- solve(t(HUk) %*% theta$Sigma2.inv %*% HUk + theta$Omega2.inv)
            beta2.mean <- beta2.variance %*% t(HUk) %*% theta$Sigma2.inv %*% data[, i]
            new_beta2[i, ] <- MASS::mvrnorm(n = 1, mu = beta2.mean, Sigma = beta2.variance)
          }
        }
        return(new_beta2)
      }
      
      proposal_function_Omega1_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        tau <- theta$tau
        n1 <- sum(tau == 2)
        Omega1.inv.df <- nu + n1 + 1
        Omega1.inv.p2 <- matrix(0, nrow = K, ncol = K)
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new.Omega1.inv <- theta$Omega1.inv
        } else{
          for (i in subj1) {
            Omega1.inv.p2 <- Omega1.inv.p2 + (theta$beta1[i, ]) %*% t(theta$beta1[i, ])
          }
          Omega1.inv.mat <- solve(Lk + Omega1.inv.p2)
          new.Omega1.inv <- MCMCpack::rwish(v = Omega1.inv.df, S = Omega1.inv.mat)
        }
        return(new.Omega1.inv)
      }
      
      proposal_function_Omega2_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        tau <- theta$tau
        n2 <- sum(tau == 1)
        Omega2.inv.df <- nu + n2 + 1
        Omega2.inv.p2 <- matrix(0, nrow = K, ncol = K)
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new.Omega2.inv <- theta$Omega2.inv
        } else{
          for (i in subj2) {
            Omega2.inv.p2 <- Omega2.inv.p2 + (theta$beta2[i, ]) %*% t(theta$beta2[i, ])
          }
          Omega2.inv.mat <- solve(Lk + Omega2.inv.p2)
          new.Omega2.inv <- MCMCpack::rwish(v = Omega2.inv.df, S = Omega2.inv.mat)
        }
        return(new.Omega2.inv)
      }
      
      proposal_function_Zi1 <- function(update.all = FALSE) {
        new.Zi1 <- theta$Zi1
        tau <- theta$tau
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new.Zi1 <- theta$Zi1
        } else{
          if (update.all){
            for (i in subj1){
              Ai1.temp <- diag(theta$D1) %*% theta$Sigma1.inv
              Ai1 <- diag(N.tps) + Ai1.temp %*% diag(theta$D1)
              ai1 <- Ai1.temp %*% (data[, i] - HUk %*% theta$beta1[i, ])
              Zi1.varcov <- solve(Ai1)
              Zi1.mean <- Zi1.varcov %*% ai1
              new.Zi1[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi1.varcov), mean = Zi1.mean, 
                                               sigma = Zi1.varcov, lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
            }
          } else{
            # Sample one value from component 1 (mvn)
            i <- sample(subj1, 1)
            Ai1.temp <- diag(theta$D1) %*% theta$Sigma1.inv
            Ai1 <- diag(N.tps) + Ai1.temp %*% diag(theta$D1)
            ai1 <- Ai1.temp %*% (data[, i] - HUk %*% theta$beta1[i, ])
            Zi1.varcov <- solve(Ai1)
            Zi1.mean <- Zi1.varcov %*% ai1
            new.Zi1[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi1.varcov), mean = Zi1.mean, 
                                             sigma = Zi1.varcov, lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
          }
        }
        
        return(new.Zi1)
      }
      
      proposal_function_Zi2 <- function(update.all = FALSE) {
        new.Zi2 <- theta$Zi2
        wi <- theta$wi
        tau <- theta$tau
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new.Zi2 <- theta$Zi2
        } else{
          if (update.all){
            for (i in subj2){
              Ai2.temp <- diag(theta$D2) %*% theta$Sigma2.inv
              Ai2 <- diag(N.tps) + Ai2.temp %*% diag(theta$D2)
              ai2 <- Ai2.temp %*% (data[, i] - HUk %*% theta$beta2[i, ])
              Zi2.varcov <- solve(Ai2)
              Zi2.mean <- Zi2.varcov %*% ai2
              new.Zi2[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi2.varcov), mean = Zi2.mean, 
                                               sigma = Zi2.varcov, lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
            }
          } else{
            # Sample one value from component 2 (mvt)
            i <- sample(subj2, 1)
            Ai2.temp <- diag(theta$D2) %*% theta$Sigma2.inv
            Ai2 <- diag(N.tps) + Ai2.temp %*% diag(theta$D2)
            ai2 <- Ai2.temp %*% (data[, i] - HUk %*% theta$beta2[i, ])
            Zi2.varcov <- solve(Ai2)
            Zi2.mean <- Zi2.varcov %*% ai2
            new.Zi2[i, ] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi2.varcov), mean = Zi2.mean, 
                                             sigma = Zi2.varcov, lower = rep(0, N.tps), upper = rep(Inf, N.tps))$samp
          }
        }
        
        return(new.Zi2)
      }
      
      proposal_function_D1 <- function() {
        tau <- theta$tau
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new.D1 <- theta$D1
        } else{
          B1.temp <- matrix(data = 0, nrow = N.tps, ncol = N.tps) + diag(x = 1/gamma, nrow = N.tps)
          b1.temp <- matrix(data = 0, nrow = N.tps, ncol = 1)
          for (i in subj1){
            B1.temp <- B1.temp + diag(theta$Zi1[i, ]) %*% theta$Sigma1.inv %*% diag(theta$Zi1[i, ])
            b1.temp <- b1.temp + diag(theta$Zi1[i, ]) %*% theta$Sigma1.inv %*%
              (data[, i] - HUk %*% theta$beta1[i,])
          }
          vecD1.mean <- solve(B1.temp) %*% b1.temp
          vecD1.varcov <- solve(B1.temp)
          new.D1 <- mnormt::rmnorm(n = 1, mean = vecD1.mean, varcov = vecD1.varcov)
        }
        return(new.D1)
      }
      
      proposal_function_D2 <- function() {
        tau <- theta$tau
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new.D2 <- theta$D2
        } else{
          B2.temp <- matrix(data = 0, nrow = N.tps, ncol = N.tps) + diag(x = 1/gamma, nrow = N.tps)
          b2.temp <- matrix(data = 0, nrow = N.tps, ncol = 1)
          for (i in subj2){
            B2.temp <- B2.temp + diag(theta$Zi2[i, ]) %*% theta$Sigma2.inv %*% diag(theta$Zi2[i, ])
            b2.temp <- b2.temp + diag(theta$Zi2[i, ]) %*% theta$Sigma2.inv %*%
              (data[, i] - HUk %*% theta$beta2[i,])
          }
          vecD2.mean <- solve(B2.temp) %*% b2.temp
          vecD2.varcov <- solve(B2.temp)
          new.D2 <- mnormt::rmnorm(n = 1, mean = vecD2.mean, varcov = vecD2.varcov)
        }
        return(new.D2)
      }
      
      proposal_function_Sigma1_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        tau <- theta$tau
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new.Sigma1.inv = theta$Sigma1.inv
        } else{
          n1 <- sum(tau == 2)
          Sigma1.inv.df <- twor + n1
          Sigma1.inv.mat <- solve(2 * kappa)
          for (i in subj1) {
            temp <- data[, i] - HUk %*% theta$beta1[i, ] - diag(theta$D1) %*% theta$Zi1[i, ]
            Sigma1.inv.mat <- Sigma1.inv.mat + temp %*% t(temp)
          }
          new.Sigma1.inv <- MCMCpack::rwish(v = Sigma1.inv.df, S = solve(Sigma1.inv.mat))
        }
        return(new.Sigma1.inv)
      }
      
      proposal_function_Sigma2_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        tau <- theta$tau
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new.Sigma2.inv = theta$Sigma2.inv
        } else{
          n2 <- sum(tau == 1)
          Sigma2.inv.df <- twor + n2
          Sigma2.inv.mat <- solve(2 * kappa)
          for (i in subj2) {
            temp <- data[, i] - HUk %*% theta$beta2[i, ] - diag(theta$D2) %*% theta$Zi2[i, ]
            Sigma2.inv.mat <- Sigma2.inv.mat + temp %*% t(temp)
          }
          new.Sigma2.inv <- MCMCpack::rwish(v = Sigma2.inv.df, S = solve(Sigma2.inv.mat))
        }
        return(new.Sigma2.inv)
      }
      
      proposal_function_wi <- function() {
        new.wi <- theta$wi
        tau <- theta$tau
        ## sample wi from a Gamma distribution
        vi <- LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1)
        wi.temp <- apply(data, 2, var)
        for (i in 1:N.subj){
          if (tau[i] == 1){ #mvt
            wi.shape <- vi[i]/2 + N.tps/2
            wi.rate <- vi[i]/2 +  N.tps * wi.temp[i]/ 2
            new.wi[i] <- rgamma(1, shape = wi.shape, rate = wi.rate)
          } else{ #mvn
            new.wi[i] <- 1
          }
        }
        return(new.wi)
      }
      
      proposal_function_pi1 <- function(){
        tau <- theta$tau
        n1 <- sum(tau == 2)
        n2 <- sum(tau == 1)
        new.pi1 <- rbeta(1, 1+n1, 1+n2)
        return(new.pi1)
      }
      
      proposal_function_tau <- function(){
        new.tau <- theta$tau
        tau <- theta$tau
        pi1 <- theta$pi1
        for (i in 1:N.subj){
          tau.A <- pi1 * LaplacesDemon::dmvn(x = data[, i],
                                             mu = t(HUk %*% theta$beta1[i,] +
                                                      diag(theta$D1) %*% theta$Zi1[i, ]),
                                             Sigma = solve(theta$Sigma1.inv))
          tau.B <- (1 - pi1) * LaplacesDemon::dmvn(x = data[, i],
                                                   mu = t(HUk %*% theta$beta2[i,] +
                                                            diag(theta$D2) %*% theta$Zi2[i, ]),
                                                   Sigma = solve(theta$Sigma2.inv)/theta$wi[i])
          if (tau.A + tau.B == 0){
            new.tau[i] <- LaplacesDemon::rbern(1, 0) + 1
          } else{
            new.tau[i] <- LaplacesDemon::rbern(1, tau.A/(tau.A + tau.B)) + 1
          }
        }
        
        new.tau <- ReMap.mm(new.tau)
        
        return(new.tau)
      }
      
      # Map each parameter to its corresponding proposal function
      proposal_functions <- list(
        beta1 = proposal_function_beta1,
        Omega1.inv = proposal_function_Omega1_inv,
        Zi1 = proposal_function_Zi1,
        D1 =  proposal_function_D1,
        Sigma1.inv = proposal_function_Sigma1_inv,
        beta2 = proposal_function_beta2,
        Omega2.inv = proposal_function_Omega2_inv,
        Zi2 = proposal_function_Zi2,
        D2 =  proposal_function_D2,
        Sigma2.inv = proposal_function_Sigma2_inv,
        wi = proposal_function_wi,
        pi1 = proposal_function_pi1,
        tau = proposal_function_tau
      )
      
      # Update the chosen parameter using its respective proposal function
      if (chosen_param == "tau"){
        theta_new[["tau"]] <- proposal_functions[["tau"]]()
        theta_new[["beta1"]] <- proposal_functions[["beta1"]](update.all = TRUE)
        theta_new[["beta2"]] <- proposal_functions[["beta2"]](update.all = TRUE)
        theta_new[["Omega1.inv"]] <- proposal_functions[["Omega1.inv"]]()
        theta_new[["Omega2.inv"]] <- proposal_functions[["Omega2.inv"]]()
        theta_new[["Zi1"]] <- proposal_functions[["Zi1"]](update.all = TRUE)
        theta_new[["Zi2"]] <- proposal_functions[["Zi2"]](update.all = TRUE)
        theta_new[["D1"]] <- proposal_functions[["D1"]]()
        theta_new[["D2"]] <- proposal_functions[["D2"]]()
        theta_new[["Sigma1.inv"]] <- proposal_functions[["Sigma1.inv"]]()
        theta_new[["Sigma2.inv"]] <- proposal_functions[["Sigma2.inv"]]()
        theta_new[["wi"]] <- proposal_functions[["wi"]]()
        theta_new[["pi1"]] <- proposal_functions[["pi1"]]()
      } else{
        theta_new[[chosen_param]] <- proposal_functions[[chosen_param]]()
      }
      
      # compute numerator and denominator of the MH ratio
      num <- (likelihood(theta_new)) * alpha[[r]] +  prior(theta_new)
      den <- (likelihood(theta)) * alpha[[r]] + prior(theta)
      
      # compute ratio
      ratio <- min(1, exp(num - den))
      
      # accept/reject step
      if(runif(1) < ratio){
        theta = theta_new
      } else{
        theta = theta
      }
      
      return(theta)
    }
    
    particles <- parLapply(cl, previous_particles, Gibbs_HM_Move)
    
    # compute the ESS
    log_incremental_w <- alphaDiff * u
    logw <- log_incremental_w + logw  # log unnormalized weights
    logmax <- max(logw)
    logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  
    W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
    logW <- log(W[[r]])                                # log normalized weights
    cat("logZ:",logZ[[r]],"\n")
    ESS[[r]] <- rESS(logW)
    
    # resample if ESS below threshold
    if( ESS[[r]]<tuning_param$eps )
    {
      cat("Resample: ESS=", ESS[[r]], '\n')
      ancestors <- systematic_resample( W[[r]] )
      particles <-  particles[ancestors]
      W[[r]] <- rep(1/N.particles,N.particles)
      logW <- log(W[[r]])
      w <- rep(1,N.particles)
      logw <- rep(0,N.particles)
    }
    previous_particles <- particles
    
  }
  
  stopCluster(cl)
  
  return(list(particles = particles,
              alpha = alpha,
              ESS = ESS,
              logZ = logZ,
              W = W,
              H = H, Uk = Uk, Lk = Lk))
}


RBFPCA.ASMC.mvn.sparse <- function(data.table, data, num.basis, num.PC, 
                                   tuning_param, seed, cov.index, search.bw,
                                   N.particles, n.cores) {
  # data.table: Input data in a matrix form
  # data: Input data in a list form
  # num.basis: Number of basis functions (J)
  # num.PC: Number of principal components (K)
  # tuning_param: List of tuning parameters for the ASMC algorithm: 1) eps - resampling threshold,
  #  2) phi - target rCESS for adaptive resampling scheme, 3) alhpa - annealing scheme. Only supply one of phi or alpha, not both.
  # seed: Seed for random number generation
  # cov.index: Covariance function index
  # search.bw: Search bandwidth
  # N.particles: Number of particles for the SMC algorithm
  # n.cores: Number of cores for the SMC algorithm
  
  set.seed(seed)
  
  # Get inputs -----------------------------
  N.subj <- length(data$x.demean)  # Number of subjects
  
  J <- num.basis  # Number of basis functions
  K <- num.PC     # Number of principal components
  
  # Define the hyperparameters -----------------------------
  CovFun <- function(s, t, cov.index){
    if (cov.index == 1){
      exp(-3*(t-s)^2)        # Exponential covariance function
    } else if (cov.index == 2){
      min(s + 1, t + 1)      # Linear covariance function
    } else{
      exp(-1*(t-s)^2)        # Exponential covariance function
    }
  }
  
  get.corner <- function(x){
    c.end <- ncol(x)
    r.end <- nrow(x)
    a = x[1,1]
    b = x[1,c.end]
    c = x[r.end,1]
    d = x[r.end,c.end]
    matrix(c(a,b,c,d), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  prior <- function(theta) {
    
    log_prior <- log(MCMCpack::dwish(theta$Omega.inv, v = nu, S = Lk.inv))
    Omega <- solve(theta$Omega.inv)  # Precompute the inverse once
    mean_zero <- rep(0, K)  # Precompute repeated vectors
    for (i in 1:N.subj){
      log_prior <- log_prior + LaplacesDemon::dwishart(Omega = theta$Sigma.inv[[id.all[i]]], nu = twor.all[i],
                                                       S = 2 * kappa.list[[i]], log = TRUE)
      
      mu_zi <- rep(0, N.tps.all[i])
      diag_sigma <- diag(1, N.tps.all[i])
      lb <- rep(0, N.tps.all[i])
      ub <- rep(Inf, N.tps.all[i])
      log_prior <- log_prior + sum(mvtnorm::dmvnorm(theta$beta[i,], mean_zero, Omega, log = TRUE))
      
      log_prior <- log_prior + sum(TruncatedNormal::dtmvnorm(theta$Zi[[id.all[i]]], mu = mu_zi, 
                                                             sigma = diag_sigma, lb = lb, ub = ub, log = TRUE))
      
      log_prior <- log_prior + sum(dnorm(theta$D[[id.all[i]]], mean = 0, sd = gamma, log = TRUE))
    }
    return(log_prior)
  }
  
  likelihood <- function(theta){
    log_likelihood <- 0
    for (i in 1:N.subj) {
      mean_part <- HUk.list[[i]] %*% theta$beta[i, ]
      zi <- theta$Zi[[id.all[i]]]
      d <- diag(theta$D[[id.all[i]]])
      sigma_inv <- theta$Sigma.inv[[id.all[i]]]
      Sigma <- solve(sigma_inv)
      log_likelihood <- log_likelihood + mvtnorm::dmvnorm(data$x.demean[[i]], 
                                                          mean = mean_part + d %*% zi, 
                                                          sigma = Sigma, log = TRUE)
    }
    return(log_likelihood)
  }
  
  # Define initial particles -----------------------------
  
  # initialize storage list
  alpha <- list() # a list for tempering parameters
  ESS <- list()   # a list for ESS
  logZ <- list()  # list for log-normalizing constant
  W <- list()
  
  # check if adaptive annealing scheme is being used
  adaptive_alpha <- "phi" %in% names(tuning_param)
  
  # initialize values
  r <- 1                  # SMC iteration
  logZ[[r]] <- 0
  
  if(adaptive_alpha){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(tuning_param$alpha)
  }
  
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl,{
    library(MCMCpack)
    library(tmvnsim)
    library(mvtnorm)
    library(mnormt)
    library(polynom)
    library(MASS)
    library(TruncatedNormal)
    library(LaplacesDemon)
  })
  
  H.list <- list()
  Uk.list <- list()
  Lk.list <- list()
  HUk.list <- list()
  kappa.list <- list()
  N.tps.all = id.all = twor.all <- c()
  
  ## nu ##
  nu <- 2*K
  
  ## a, b ##
  a = b <- 10^(-3)
  
  ## theta_0 ##
  theta0 <- rep(0, K)
  
  ## gamma ##
  gamma <- 10 
  
  for (i in 1:N.subj){
    
    N.tps <- length(data$pp[[i]]) # number of observations for this subject
    N.tps.all[i] <- N.tps
    this.id <- names(data$pp)[[i]]  # id for this subject
    id.all[i] <- this.id
    
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
    HUk = H %*% Uk
    
    H.list[[i]] <- H
    Uk.list[[i]] <- Uk
    HUk.list[[i]] <- HUk
    Lk.list[[i]] <- Lk
    
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
    twor.all[i] <- N.tps
    
    ## kappa ##
    kappa <- 100 / twor * inv(R)
    kappa.list[[i]] <- kappa
    
  }
  
  Lk.inv <- solve(apply(simplify2array(Lk.list), 1:2, mean))
  
  # initialize particles
  init <- function(j)
  {
    set.seed(j)
    init_particle <- list(
      beta = matrix(0, nrow = N.subj, ncol = K),
      Omega.inv = MCMCpack::rwish(v = nu, S = Lk.inv), 
      Zi = list(),
      D = list(),
      Sigma.inv = list()
    )
    
    is_invertible <- function(X) !inherits(try(solve(X), silent = TRUE), "try-error")
    
    init_particle$Omega.inv <- MCMCpack::rwish(v = twor.all[i], S = Lk.inv)
    while (!is_invertible(init_particle$Omega.inv)){
      init_particle$Omega.inv <- MCMCpack::rwish(v = twor.all[i], S = Lk.inv)
    }
    Omega <- solve(init_particle$Omega.inv)
    
    for (i in 1:N.subj) {
      init_particle$beta[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega)
      
      init_particle$Zi[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps.all[i], 
                                                        mean = rep(0,N.tps.all[i]), 
                                                        sigma = diag(1,N.tps.all[i]), 
                                                        lower = rep(0, N.tps.all[i]), 
                                                        upper = rep(Inf, N.tps.all[i]))$samp
      init_particle$D[[id.all[i]]] <- rnorm(n = N.tps.all[i], mean = 0, sd = gamma)
      init_particle$Sigma.inv[[id.all[i]]] <- MCMCpack::rwish(v = N.tps.all[i], S = 2*kappa.list[[i]])
    }
    return(init_particle)
  }

  previous_particles <- parLapply(cl, 1:N.particles, function(k) init(k))
  
  ## Initialize weights
  # normalized weights
  # Particles were sampled from the same reference distribution. They have the equal weights.
  W[[r]] <- rep(1/N.particles,N.particles)
  logW <- log(W[[r]])
  
  # unnormalized weights
  w <- rep(1,N.particles)
  logw <- rep(0,N.particles)  
  
  # Annealed SMC Algorithm -----------------------------
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    # evaluate the log-likelihood for updating alpha
    u <- rep(0, N.particles)   # incremental log-importance weights
    u <- sapply(1:N.particles, function(k){
      logL <- likelihood(previous_particles[[k]])
      return(logL)
    })
    
    if(adaptive_alpha){
      # update alpha with bisection
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$phi )
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n")
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]]
    }
    
    Gibbs_HM_Move <- function(theta)
    {
      theta_new <- theta
      # Randomly pick one parameter to update
      chosen_param <- sample(names(theta_new), 1)
      
      proposal_function_beta <- function() {
        new_beta <- theta$beta
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        beta.variance <- solve(t(HUk.list[[i]]) %*% theta$Sigma.inv[[id.all[i]]] %*% HUk.list[[i]] + theta$Omega.inv)
        beta.mean <- beta.variance %*% t(HUk.list[[i]]) %*% theta$Sigma.inv[[id.all[i]]] %*% data$x.demean[[i]]
        new_beta[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
        return(new_beta)
      }
      
      proposal_function_Omega_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        new.Omega.inv <- theta$Omega.inv
        Omega.inv.df <- nu + N.subj + 1
        Omega.inv.p2 <- matrix(0, nrow = K, ncol = K)
        for (i in 1:N.subj) {
          Omega.inv.p2 <- Omega.inv.p2 + (theta$beta[i, ]) %*% t(theta$beta[i, ])
        }
        Omega.inv.mat <- solve(Lk + Omega.inv.p2)
        new.Omega.inv <- MCMCpack::rwish(v = Omega.inv.df, S = Omega.inv.mat)
        return(new.Omega.inv)
      }
      
      proposal_function_Zi <- function() {
        new.Zi <- theta$Zi
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        Ai.temp <- diag(theta$D[[id.all[i]]]) %*% theta$Sigma.inv[[id.all[i]]]
        Ai <- diag(N.tps.all[i]) + Ai.temp %*% diag(theta$D[[id.all[i]]])
        ai <- Ai.temp %*% (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta[i, ])
        Zi.varcov <- solve(Ai)
        Zi.mean <- Zi.varcov %*% ai
        new.Zi[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), 
                                                mean = Zi.mean, sigma = Zi.varcov, 
                                                lower = rep(0, N.tps.all[i]), 
                                                upper = rep(Inf, N.tps.all[i]))$samp
        return(new.Zi)
      }
      
      proposal_function_D <- function() {
        new.D <- theta$D
        for (i in 1:N.subj){
          B <- diag(x = 1/gamma, nrow = N.tps.all[i]) + diag(theta$Zi[[id.all[i]]]) %*% 
            theta$Sigma.inv[[id.all[i]]] %*% diag(theta$Zi[[id.all[i]]])
          b <- diag(theta$Zi[[id.all[i]]]) %*% theta$Sigma.inv[[id.all[i]]] %*%
            (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta[i, ])
          vecD.mean <- solve(B) %*% b
          vecD.varcov <- solve(B)
          new.D[[id.all[i]]] <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
        }
        return(new.D)
      }
      
      proposal_function_Sigma_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        new.Sigma.inv <- theta$Sigma.inv
        for (i in 1:N.subj) {
          Sigma.inv.df <- twor.all[i] + N.subj
          Sigma.inv.mat <- solve(2 * kappa.list[[i]])
          temp <- data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta[i, ] - diag(theta$D[[id.all[i]]]) %*% theta$Zi[[id.all[i]]]
          Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
          new.Sigma.inv[[id.all[i]]] <- MCMCpack::rwish(v = Sigma.inv.df, S = solve(Sigma.inv.mat))
        }
        return(new.Sigma.inv)
      }
      
      # Map each parameter to its corresponding proposal function
      proposal_functions <- list(
        beta = proposal_function_beta,
        Omega.inv = proposal_function_Omega_inv,
        Zi = proposal_function_Zi,
        D =  proposal_function_D,
        Sigma.inv = proposal_function_Sigma_inv
      )
      
      # Update the chosen parameter using its respective proposal function
      theta_new[[chosen_param]] <- proposal_functions[[chosen_param]]()
      
      # compute numerator and denominator of the MH ratio
      num <- (likelihood(theta_new)) * alpha[[r]] +  prior(theta_new)
      den <- (likelihood(theta)) * alpha[[r]] + prior(theta)
      
      # compute ratio
      ratio <- min(1, exp(num - den))
      
      # accept/reject step
      if(runif(1) < ratio){
        theta = theta_new
      } else{
        theta = theta
      }
      
      return(theta)
    }
    
    particles <- parLapply(cl, previous_particles, Gibbs_HM_Move)
    
    # compute the ESS
    log_incremental_w <- alphaDiff * u
    logw <- log_incremental_w + logw  # log unnormalized weights
    logmax <- max(logw)
    logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  
    W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
    logW <- log(W[[r]])                                # log normalized weights
    cat("logZ:",logZ[[r]],"\n")
    ESS[[r]] <- rESS(logW)
    
    # resample if ESS below threshold
    if( ESS[[r]]<tuning_param$eps )
    {
      cat("Resample: ESS=", ESS[[r]], '\n')
      ancestors <- systematic_resample( W[[r]] )
      particles <-  particles[ancestors]
      W[[r]] <- rep(1/N.particles,N.particles)
      logW <- log(W[[r]])
      w <- rep(1,N.particles)
      logw <- rep(0,N.particles)
    }
    previous_particles <- particles
    
  }
  stopCluster(cl)
  
  return(list(particles = particles,
              alpha = alpha,
              ESS = ESS,
              logZ = logZ,
              W = W,
              H = H.list, Uk = Uk.list, Lk = Lk.list))
}


RBFPCA.ASMC.mvt.sparse <- function(data.table, data, num.basis, num.PC, 
                                   tuning_param, seed, cov.index, search.bw,
                                   N.particles, n.cores) {
  # data.table: Input data in a matrix form
  # data: Input data in a list form
  # num.basis: Number of basis functions (J)
  # num.PC: Number of principal components (K)
  # tuning_param: List of tuning parameters for the ASMC algorithm: 1) eps - resampling threshold,
  #  2) phi - target rCESS for adaptive resampling scheme, 3) alhpa - annealing scheme. Only supply one of phi or alpha, not both.
  # seed: Seed for random number generation
  # cov.index: Covariance function index
  # search.bw: Search bandwidth
  # N.particles: Number of particles for the SMC algorithm
  # n.cores: Number of cores for the SMC algorithm
  
  set.seed(seed)
  
  # Get inputs -----------------------------
  N.subj <- length(data$x.demean)  # Number of subjects
  
  J <- num.basis  # Number of basis functions
  K <- num.PC     # Number of principal components
  
  # Define the hyperparameters -----------------------------
  CovFun <- function(s, t, cov.index){
    if (cov.index == 1){
      exp(-3*(t-s)^2)        # Exponential covariance function
    } else if (cov.index == 2){
      min(s + 1, t + 1)      # Linear covariance function
    } else{
      exp(-1*(t-s)^2)        # Exponential covariance function
    }
  }
  
  get.corner <- function(x){
    c.end <- ncol(x)
    r.end <- nrow(x)
    a = x[1,1]
    b = x[1,c.end]
    c = x[r.end,1]
    d = x[r.end,c.end]
    matrix(c(a,b,c,d), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  prior <- function(theta) {
    
    log_prior <- log(MCMCpack::dwish(theta$Omega.inv, v = nu, S = Lk.inv))
    log_prior <- log_prior + sum(LaplacesDemon::dtrunc(theta$v, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1))
    Omega <- solve(theta$Omega.inv)  # Precompute the inverse once
    mean_zero <- rep(0, K)  # Precompute repeated vectors
    for (i in 1:N.subj){
      log_prior <- log_prior + LaplacesDemon::dwishart(Omega = theta$Sigma.inv[[id.all[i]]], nu = twor.all[i],
                                                       S = 2 * kappa.list[[i]], log = TRUE)
      
      mu_zi <- rep(0, N.tps.all[i])
      diag_sigma <- diag(1, N.tps.all[i])
      lb <- rep(0, N.tps.all[i])
      ub <- rep(Inf, N.tps.all[i])
      log_prior <- log_prior + sum(mvtnorm::dmvnorm(theta$beta[i,], mean_zero, Omega, log = TRUE))
      
      log_prior <- log_prior + sum(TruncatedNormal::dtmvnorm(theta$Zi[[id.all[i]]], mu = mu_zi, 
                                                             sigma = diag_sigma, lb = lb, ub = ub, log = TRUE))
      
      log_prior <- log_prior + sum(dnorm(theta$D[[id.all[i]]], mean = 0, sd = gamma, log = TRUE))
      
      # Sum wi for MVT
      log_prior <- log_prior + sum(dgamma(theta$wi[i], shape = theta$v/2, rate = theta$v/2))
      
    }
    return(log_prior)
  }
  
  likelihood <- function(theta){
    log_likelihood <- 0
    for (i in 1:N.subj) {
      mean_part <- HUk.list[[i]] %*% theta$beta[i, ]
      zi <- theta$Zi[[id.all[i]]]
      d <- diag(theta$D[[id.all[i]]])
      sigma_inv <- theta$Sigma.inv[[id.all[i]]]
      Sigma <- solve(sigma_inv)
      wi <- theta$wi[i]
      log_likelihood <- log_likelihood + mvtnorm::dmvnorm(data$x.demean[[i]], 
                                                          mean = mean_part + d %*% zi, 
                                                          sigma = Sigma/wi, log = TRUE)
    }
    return(log_likelihood)
  }
  
  # Define initial particles -----------------------------
  
  # initialize storage list
  alpha <- list() # a list for tempering parameters
  ESS <- list()   # a list for ESS
  logZ <- list()  # list for log-normalizing constant
  W <- list()
  
  # check if adaptive annealing scheme is being used
  adaptive_alpha <- "phi" %in% names(tuning_param)
  
  # initialize values
  r <- 1                  # SMC iteration
  logZ[[r]] <- 0
  
  if(adaptive_alpha){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(tuning_param$alpha)
  }
  
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl,{
    library(MCMCpack)
    library(tmvnsim)
    library(mvtnorm)
    library(mnormt)
    library(polynom)
    library(MASS)
    library(TruncatedNormal)
    library(LaplacesDemon)
  })
  
  H.list <- list()
  Uk.list <- list()
  Lk.list <- list()
  HUk.list <- list()
  kappa.list <- list()
  N.tps.all = id.all = twor.all <- c()
  
  ## nu ##
  nu <- 2*K
  
  ## a, b ##
  a = b <- 10^(-3)
  
  ## theta_0 ##
  theta0 <- rep(0, K)
  
  ## gamma ##
  gamma <- 10 
  
  for (i in 1:N.subj){
    
    N.tps <- length(data$pp[[i]]) # number of observations for this subject
    N.tps.all[i] <- N.tps
    this.id <- names(data$pp)[[i]]  # id for this subject
    id.all[i] <- this.id
    
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
    #Lk.inv <- solve(Lk)
    HUk = H %*% Uk
    
    H.list[[i]] <- H
    Uk.list[[i]] <- Uk
    HUk.list[[i]] <- HUk
    Lk.list[[i]] <- Lk
    
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
    twor.all[i] <- N.tps
    
    ## kappa ##
    kappa <- 100 / twor * inv(R)
    kappa.list[[i]] <- kappa
    
  }
  
  Lk.inv <- solve(apply(simplify2array(Lk.list), 1:2, mean))
  
  # initialize particles
  init <- function(j)
  {
    set.seed(j)
    wi_func <- function(x){
      rgamma(1,shape = x/2, rate = x/2)
    }
    init_particle <- list(
      beta = matrix(0, nrow = N.subj, ncol = K),
      Omega.inv = MCMCpack::rwish(v = nu, S = Lk.inv), 
      Zi = list(),
      D = list(),
      Sigma.inv = list(),
      vi = LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1),
      wi = c()
    )

    is_invertible <- function(X) !inherits(try(solve(X), silent = TRUE), "try-error")
    
    init_particle$Omega.inv <- MCMCpack::rwish(v = nu, S = Lk.inv)
    while (!is_invertible(init_particle$Omega.inv)){
      init_particle$Omega.inv <- MCMCpack::rwish(v = nu, S = Lk.inv)
    }
    Omega <- solve(init_particle$Omega.inv)
    
    init_particle$wi <- sapply(init_particle$vi, wi_func)
    
    for (i in 1:N.subj) {
      init_particle$beta[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega)
      
      init_particle$Zi[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps.all[i], 
                                                        mean = rep(0,N.tps.all[i]), 
                                                        sigma = diag(1,N.tps.all[i]), 
                                                        lower = rep(0, N.tps.all[i]), 
                                                        upper = rep(Inf, N.tps.all[i]))$samp
      init_particle$D[[id.all[i]]] <- rnorm(n = N.tps.all[i], mean = 0, sd = gamma)
      init_particle$Sigma.inv[[id.all[i]]] <- MCMCpack::rwish(v = twor.all[i], S = 2*kappa.list[[i]])
    }
    return(init_particle)
  }
 
  previous_particles <- parLapply(cl, 1:N.particles, function(k) init(k))
  
  ## Initialize weights
  # normalized weights
  # Particles were sampled from the same reference distribution. They have the equal weights.
  W[[r]] <- rep(1/N.particles,N.particles)
  logW <- log(W[[r]])
  
  # unnormalized weights
  w <- rep(1,N.particles)
  logw <- rep(0,N.particles)  
  
  # Annealed SMC Algorithm -----------------------------
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    # evaluate the log-likelihood for updating alpha
    u <- rep(0, N.particles)   # incremental log-importance weights
    u <- sapply(1:N.particles, function(k){
      logL <- likelihood(previous_particles[[k]])
      return(logL)
    })
    
    if(adaptive_alpha){
      # update alpha with bisection
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$phi )
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n")
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]]
    }
    
    Gibbs_HM_Move <- function(theta)
    {
      theta_new <- theta
      # Randomly pick one parameter to update
      param_names <- names(theta_new)[!names(theta_new) == "vi"]
      chosen_param <- sample(param_names, 1)
      
      proposal_function_beta <- function() {
        new_beta <- theta$beta
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        beta.variance <- solve(t(HUk.list[[i]]) %*% theta$Sigma.inv[[id.all[i]]] %*% HUk.list[[i]] + theta$Omega.inv)
        beta.mean <- beta.variance %*% t(HUk.list[[i]]) %*% theta$Sigma.inv[[id.all[i]]] %*% data$x.demean[[i]]
        new_beta[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
        return(new_beta)
      }
      
      proposal_function_Omega_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        new.Omega.inv <- theta$Omega.inv
        Omega.inv.df <- nu + N.subj + 1
        Omega.inv.p2 <- matrix(0, nrow = K, ncol = K)
        for (i in 1:N.subj) {
          Omega.inv.p2 <- Omega.inv.p2 + (theta$beta[i, ]) %*% t(theta$beta[i, ])
        }
        Omega.inv.mat <- solve(Lk + Omega.inv.p2)
        new.Omega.inv <- MCMCpack::rwish(v = Omega.inv.df, S = Omega.inv.mat)
        return(new.Omega.inv)
      }
      
      proposal_function_Zi <- function() {
        new.Zi <- theta$Zi
        # Sample one value from 1 to N.subj
        i <- sample(1:N.subj, 1)
        Ai.temp <- diag(theta$D[[id.all[i]]]) %*% theta$Sigma.inv[[id.all[i]]]
        Ai <- diag(N.tps.all[i]) + Ai.temp %*% diag(theta$D[[id.all[i]]])
        ai <- Ai.temp %*% (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta[i, ])
        Zi.varcov <- solve(Ai)
        Zi.mean <- Zi.varcov %*% ai
        new.Zi[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), 
                                                mean = Zi.mean, sigma = Zi.varcov, 
                                                lower = rep(0, N.tps.all[i]), 
                                                upper = rep(Inf, N.tps.all[i]))$samp
        return(new.Zi)
      }
      
      proposal_function_D <- function() {
        new.D <- theta$D
        for (i in 1:N.subj){
          B <- diag(x = 1/gamma, nrow = N.tps.all[i]) + diag(theta$Zi[[id.all[i]]]) %*% 
            theta$Sigma.inv[[id.all[i]]] %*% diag(theta$Zi[[id.all[i]]])
          b <- diag(theta$Zi[[id.all[i]]]) %*% theta$Sigma.inv[[id.all[i]]] %*%
            (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta[i, ])
          vecD.mean <- solve(B) %*% b
          vecD.varcov <- solve(B)
          new.D[[id.all[i]]] <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
        }
        return(new.D)
      }
      
      proposal_function_Sigma_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        new.Sigma.inv <- theta$Sigma.inv
        for (i in 1:N.subj) {
          Sigma.inv.df <- twor.all[i] + N.subj
          Sigma.inv.mat <- solve(2 * kappa.list[[i]])
          temp <- data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta[i, ] - diag(theta$D[[id.all[i]]]) %*% theta$Zi[[id.all[i]]]
          Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
          new.Sigma.inv[[id.all[i]]] <- MCMCpack::rwish(v = Sigma.inv.df, S = solve(Sigma.inv.mat))
        }
        return(new.Sigma.inv)
      }
      
      proposal_function_wi <- function() {
        new.wi <- theta$wi
        ## sample wi from a Gamma distribution
        vi <- LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1)
        for (i in 1:N.subj){
          wi.temp <- var(data$x.demean[[i]])
          wi.shape <- vi[i]/2 + N.subj/2
          wi.rate <- vi[i]/2 +  N.subj * wi.temp/ 2
          new.wi[i] <- rgamma(1, shape = wi.shape, rate = wi.rate)
        }
        return(new.wi)
      }
      
      # Map each parameter to its corresponding proposal function
      proposal_functions <- list(
        beta = proposal_function_beta,
        Omega.inv = proposal_function_Omega_inv,
        Zi = proposal_function_Zi,
        D =  proposal_function_D,
        Sigma.inv = proposal_function_Sigma_inv,
        wi = proposal_function_wi
      )
      
      # Update the chosen parameter using its respective proposal function
      theta_new[[chosen_param]] <- proposal_functions[[chosen_param]]()
      
      # compute numerator and denominator of the MH ratio
      num <- (likelihood(theta_new)) * alpha[[r]] +  prior(theta_new)
      den <- (likelihood(theta)) * alpha[[r]] + prior(theta)
      
      # compute ratio
      ratio <- min(1, exp(num - den))
      
      # accept/reject step
      if(runif(1) < ratio){
        theta = theta_new
      } else{
        theta = theta
      }
      
      return(theta)
    }
    
    particles <- parLapply(cl, previous_particles, Gibbs_HM_Move)
    
    # compute the ESS
    log_incremental_w <- alphaDiff * u
    logw <- log_incremental_w + logw  # log unnormalized weights
    logmax <- max(logw)
    logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  
    W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
    logW <- log(W[[r]])                                # log normalized weights
    cat("logZ:",logZ[[r]],"\n")
    ESS[[r]] <- rESS(logW)
    
    # resample if ESS below threshold
    if( ESS[[r]]<tuning_param$eps )
    {
      cat("Resample: ESS=", ESS[[r]], '\n')
      ancestors <- systematic_resample( W[[r]] )
      particles <-  particles[ancestors]
      W[[r]] <- rep(1/N.particles,N.particles)
      logW <- log(W[[r]])
      w <- rep(1,N.particles)
      logw <- rep(0,N.particles)
    }
    previous_particles <- particles
    
  }
  stopCluster(cl)
  
  return(list(particles = particles,
              alpha = alpha,
              ESS = ESS,
              logZ = logZ,
              W = W,
              H = H.list, Uk = Uk.list, Lk = Lk.list))
}


RBFPCA.ASMC.mm.sparse <- function(data.table, data, num.basis, num.PC, 
                                  tuning_param, seed, cov.index, search.bw,
                                  N.particles, n.cores) {
  # data.table: Input data in a matrix form
  # data: Input data in a list form
  # num.basis: Number of basis functions (J)
  # num.PC: Number of principal components (K)
  # tuning_param: List of tuning parameters for the ASMC algorithm: 1) eps - resampling threshold,
  #  2) phi - target rCESS for adaptive resampling scheme, 3) alhpa - annealing scheme. Only supply one of phi or alpha, not both.
  # seed: Seed for random number generation
  # cov.index: Covariance function index
  # search.bw: Search bandwidth
  # N.particles: Number of particles for the SMC algorithm
  # n.cores: Number of cores for the SMC algorithm
  
  set.seed(seed)
  
  # Get inputs -----------------------------
  N.subj <- length(data$x.demean)  # Number of subjects
  
  J <- num.basis  # Number of basis functions
  K <- num.PC     # Number of principal components
  
  # Define the hyperparameters -----------------------------
  CovFun <- function(s, t, cov.index){
    if (cov.index == 1){
      exp(-3*(t-s)^2)        # Exponential covariance function
    } else if (cov.index == 2){
      min(s + 1, t + 1)      # Linear covariance function
    } else{
      exp(-1*(t-s)^2)        # Exponential covariance function
    }
  }
  
  get.corner <- function(x){
    c.end <- ncol(x)
    r.end <- nrow(x)
    a = x[1,1]
    b = x[1,c.end]
    c = x[r.end,1]
    d = x[r.end,c.end]
    matrix(c(a,b,c,d), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  prior <- function(theta) {
    
    log_prior <- log(MCMCpack::dwish(theta$Omega1.inv, v = nu, S = Lk.inv)) +
      log(MCMCpack::dwish(theta$Omega1.inv, v = nu, S = Lk.inv))
    log_prior <- log_prior + sum(LaplacesDemon::dtrunc(theta$v, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1))
    Omega1 <- solve(theta$Omega1.inv)  # Precompute the inverse once
    Omega2 <- solve(theta$Omega2.inv)  # Precompute the inverse once
    mean_zero <- rep(0, K)  # Precompute repeated vectors
    tau <- theta$tau
    for (i in 1:N.subj){
      if (tau[i] == 1){
        log_prior <- log_prior + LaplacesDemon::dwishart(Omega = theta$Sigma2.inv[[id.all[i]]], nu = twor.all[i],
                                                         S = 2 * kappa.list[[i]], log = TRUE)
        mu_zi <- rep(0, N.tps.all[i])
        diag_sigma <- diag(1, N.tps.all[i])
        lb <- rep(0, N.tps.all[i])
        ub <- rep(Inf, N.tps.all[i])
        log_prior <- log_prior + sum(mvtnorm::dmvnorm(theta$beta2[i,], mean_zero, Omega2, log = TRUE))
        
        log_prior <- log_prior + sum(TruncatedNormal::dtmvnorm(theta$Zi2[[id.all[i]]], mu = mu_zi, 
                                                               sigma = diag_sigma, lb = lb, ub = ub, log = TRUE))
        
        log_prior <- log_prior + sum(dnorm(theta$D2[[id.all[i]]], mean = 0, sd = gamma, log = TRUE))
        
        # Sum wi for MVT
        log_prior <- log_prior + sum(dgamma(theta$wi[i], shape = theta$v[i]/2, rate = theta$v[i]/2))
      } else{
        log_prior <- log_prior + LaplacesDemon::dwishart(Omega = theta$Sigma1.inv[[id.all[i]]], nu = twor.all[i],
                                                         S = 2 * kappa.list[[i]], log = TRUE)
        mu_zi <- rep(0, N.tps.all[i])
        diag_sigma <- diag(1, N.tps.all[i])
        lb <- rep(0, N.tps.all[i])
        ub <- rep(Inf, N.tps.all[i])
        log_prior <- log_prior + sum(mvtnorm::dmvnorm(theta$beta1[i,], mean_zero, Omega1, log = TRUE))
        
        log_prior <- log_prior + sum(TruncatedNormal::dtmvnorm(theta$Zi1[[id.all[i]]], mu = mu_zi, 
                                                               sigma = diag_sigma, lb = lb, ub = ub, log = TRUE))
        
        log_prior <- log_prior + sum(dnorm(theta$D1[[id.all[i]]], mean = 0, sd = gamma, log = TRUE))
      }
    }
    return(log_prior)
  }
  
  likelihood <- function(theta){
    log_likelihood <- 0
    tau <-  theta$tau
    pi1 <- theta$pi1

    for (i in 1:N.subj) {
      if (tau[i] == 1){
        mean_part <- HUk.list[[i]] %*% theta$beta2[i, ]
        zi <- theta$Zi2[[id.all[i]]]
        d <- diag(theta$D2[[id.all[i]]])
        sigma_inv <- theta$Sigma2.inv[[id.all[i]]]
        Sigma <- solve(sigma_inv)
        wi <- theta$wi[i]
        log_likelihood <- log_likelihood + log(1-pi1) +  mvtnorm::dmvnorm(data$x.demean[[i]], 
                                                                          mean = mean_part + d %*% zi, 
                                                                          sigma = Sigma/wi, log = TRUE)
      } else{
        mean_part <- HUk.list[[i]] %*% theta$beta1[i, ]
        zi <- theta$Zi1[[id.all[i]]]
        d <- diag(theta$D1[[id.all[i]]])
        sigma_inv <- theta$Sigma1.inv[[id.all[i]]]
        Sigma <- solve(sigma_inv)
        log_likelihood <- log_likelihood + log(pi1) +  mvtnorm::dmvnorm(data$x.demean[[i]], 
                                                                        mean = mean_part + d %*% zi, 
                                                                        sigma = Sigma, log = TRUE)
      }
    }
    return(log_likelihood)
  }
  
  # Define initial particles -----------------------------
  
  # initialize storage list
  alpha <- list() # a list for tempering parameters
  ESS <- list()   # a list for ESS
  logZ <- list()  # list for log-normalizing constant
  W <- list()
  
  # check if adaptive annealing scheme is being used
  adaptive_alpha <- "phi" %in% names(tuning_param)
  
  # initialize values
  r <- 1                  # SMC iteration
  logZ[[r]] <- 0
  
  if(adaptive_alpha){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(tuning_param$alpha)
  }
  
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl,{
    library(MCMCpack)
    library(tmvnsim)
    library(mvtnorm)
    library(mnormt)
    library(polynom)
    library(MASS)
    library(TruncatedNormal)
    library(LaplacesDemon)
  })
  
  H.list <- list()
  Uk.list <- list()
  Lk.list <- list()
  HUk.list <- list()
  kappa.list <- list()
  N.tps.all = id.all = twor.all <- c()
  
  ## nu ##
  nu <- 2*K
  
  ## a, b ##
  a = b <- 10^(-3)
  
  ## theta_0 ##
  theta0 <- rep(0, K)
  
  ## gamma ##
  gamma <- 10 
  
  for (i in 1:N.subj){
    
    N.tps <- length(data$pp[[i]]) # number of observations for this subject
    N.tps.all[i] <- N.tps
    this.id <- names(data$pp)[[i]]  # id for this subject
    id.all[i] <- this.id
    
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
    HUk = H %*% Uk
    
    H.list[[i]] <- H
    Uk.list[[i]] <- Uk
    HUk.list[[i]] <- HUk
    Lk.list[[i]] <- Lk
    
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
    twor.all[i] <- N.tps
    
    ## kappa ##
    kappa <- 100 / twor * inv(R)
    kappa.list[[i]] <- kappa
    
  }
  
  Lk.inv <- solve(apply(simplify2array(Lk.list), 1:2, mean))
  
  # initialize particles
  init <- function(j)
  {
    set.seed(j)
    ReMap.mm <- function(x){
      x.unique <- unique(x)
      x.remap <- vector(mode = "integer", length = length(x))
      for (i in 1:length(x.unique)){
        x.remap[x == x.unique[i]] <- i
      }
      return(x.remap)
    }
    wi_func <- function(x){
      rgamma(1,shape = x/2, rate = x/2)
    }
    init_particle <- list(
      beta1 = matrix(0, nrow = N.subj, ncol = K),
      beta2 = matrix(0, nrow = N.subj, ncol = K),
      Omega1.inv = MCMCpack::rwish(v = nu, S = Lk.inv), 
      Omega2.inv = MCMCpack::rwish(v = nu, S = Lk.inv), 
      Zi1 = list(),
      Zi2 = list(),
      D1 = list(),
      D2 = list(),
      Sigma1.inv = list(),
      Sigma2.inv = list(),
      vi = LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1),
      wi = c(),
      pi1 = rbeta(1,1,1)
    )
    
    is_invertible <- function(X) !inherits(try(solve(X), silent = TRUE), "try-error")
    
    init_particle$tau <- ReMap.mm((LaplacesDemon::rbern(N.subj, init_particle$pi1) + 1))
    while (length(unique(init_particle$tau)) == 1){
      init_particle$pi1 <- rbeta(1,1,1)
      init_particle$tau <- ReMap.mm((LaplacesDemon::rbern(N.subj, init_particle$pi1) + 1))
    }
    
    while (!is_invertible(init_particle$Omega1.inv)){
      init_particle$Omega1.inv <- MCMCpack::rwish(v = nu, S = Lk.inv)
    }
    Omega1 <- solve(init_particle$Omega1.inv)
    
    while (!is_invertible(init_particle$Omega2.inv)){
      init_particle$Omega2.inv <- MCMCpack::rwish(v = nu, S = Lk.inv)
    }
    Omega2 <- solve(init_particle$Omega2.inv)
    
    init_particle$wi <- sapply(init_particle$vi, wi_func)
    
    for (i in 1:N.subj) {
      init_particle$beta1[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega1)
      init_particle$beta2[i, ] <- MASS::mvrnorm(n = 1, mu = rep(0,K), Sigma = Omega2)
      init_particle$Zi1[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps.all[i], 
                                                         mean = rep(0,N.tps.all[i]), 
                                                         sigma = diag(1,N.tps.all[i]), 
                                                         lower = rep(0, N.tps.all[i]), 
                                                         upper = rep(Inf, N.tps.all[i]))$samp
      init_particle$Zi2[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = N.tps.all[i], 
                                                         mean = rep(0,N.tps.all[i]), 
                                                         sigma = diag(1,N.tps.all[i]), 
                                                         lower = rep(0, N.tps.all[i]), 
                                                         upper = rep(Inf, N.tps.all[i]))$samp
      init_particle$D1[[id.all[i]]] <- rnorm(n = N.tps.all[i], mean = 0, sd = gamma)
      init_particle$D2[[id.all[i]]] <- rnorm(n = N.tps.all[i], mean = 0, sd = gamma)
      init_particle$Sigma1.inv[[id.all[i]]] <- MCMCpack::rwish(v = twor.all[i], S = 2*kappa.list[[i]])
      init_particle$Sigma2.inv[[id.all[i]]] <- MCMCpack::rwish(v = twor.all[i], S = 2*kappa.list[[i]])
    }
    return(init_particle)
  }
  
  previous_particles <- parLapply(cl, 1:N.particles, function(k) init(k))
  
  ## Initialize weights
  # normalized weights
  # Particles were sampled from the same reference distribution. They have the equal weights.
  W[[r]] <- rep(1/N.particles,N.particles)
  logW <- log(W[[r]])
  
  # unnormalized weights
  w <- rep(1,N.particles)
  logw <- rep(0,N.particles)  
  
  # Annealed SMC Algorithm -----------------------------
  while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
  {
    cat("iteration:",r,"\n")
    r <- r+1  # increment iteration
    
    # evaluate the log-likelihood for updating alpha
    u <- rep(0, N.particles)   # incremental log-importance weights
    u <- sapply(1:N.particles, function(k){
      logL <- likelihood(previous_particles[[k]])
      return(logL)
    })
    
    if(adaptive_alpha){
      # update alpha with bisection
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$phi )
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n")
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]]
    }
    
    Gibbs_HM_Move <- function(theta)
    {
      
      ReMap.mm <- function(x){
        x.unique <- unique(x)
        x.remap <- vector(mode = "integer", length = length(x))
        for (i in 1:length(x.unique)){
          x.remap[x == x.unique[i]] <- i
        }
        return(x.remap)
      }
      
      theta_new <- theta
      # Randomly pick one parameter to update
      param_names <- names(theta_new)[!names(theta_new) == "vi"]
      chosen_param <- sample(param_names, 1)
      print(chosen_param)
      
      proposal_function_beta1 <- function(update.all = FALSE) {
        new_beta1 <- theta$beta1
        tau <- theta$tau
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new_beta1 <- theta$beta1
        } else{
          if (update.all){
            for (i in subj1){
              beta.variance <- solve(t(HUk.list[[i]]) %*% theta$Sigma1.inv[[id.all[i]]] %*% HUk.list[[i]] + theta$Omega1.inv)
              beta.mean <- beta.variance %*% t(HUk.list[[i]]) %*% theta$Sigma1.inv[[id.all[i]]] %*% data$x.demean[[i]]
              new_beta1[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
            }
          } else{
            # Sample one value from subj1
            i <- sample(subj1, 1)
            beta.variance <- solve(t(HUk.list[[i]]) %*% theta$Sigma1.inv[[id.all[i]]] %*% HUk.list[[i]] + theta$Omega1.inv)
            beta.mean <- beta.variance %*% t(HUk.list[[i]]) %*% theta$Sigma1.inv[[id.all[i]]] %*% data$x.demean[[i]]
            new_beta1[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
          }
        }
        return(new_beta1)
      }
      
      proposal_function_beta2 <- function(update.all = FALSE) {
        new_beta2 <- theta$beta2
        tau <- theta$tau
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new_beta2 <- theta$beta2
        } else{
          if (update.all){
            for (i in subj2){
              beta.variance <- solve(t(HUk.list[[i]]) %*% theta$Sigma2.inv[[id.all[i]]] %*% HUk.list[[i]] + theta$Omega2.inv)
              beta.mean <- beta.variance %*% t(HUk.list[[i]]) %*% theta$Sigma2.inv[[id.all[i]]] %*% data$x.demean[[i]]
              new_beta2[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
            }
          } else{
            # Sample one value from subj2
            i <- sample(subj2, 1)
            beta.variance <- solve(t(HUk.list[[i]]) %*% theta$Sigma2.inv[[id.all[i]]] %*% HUk.list[[i]] + theta$Omega2.inv)
            beta.mean <- beta.variance %*% t(HUk.list[[i]]) %*% theta$Sigma2.inv[[id.all[i]]] %*% data$x.demean[[i]]
            new_beta2[i, ] <- MASS::mvrnorm(n = 1, mu = beta.mean, Sigma = beta.variance)
          }
        }
        return(new_beta2)
      }
      
      proposal_function_Omega1_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        tau <- theta$tau
        n1 <- sum(tau == 2)
        subj1 <- which(tau == 2)
        new.Omega1.inv <- theta$Omega1.inv
        Omega.inv.df <- nu + n1 + 1
        Omega.inv.p2 <- matrix(0, nrow = K, ncol = K)
        if (length(subj1) == 0){
          new.Omega1.inv <- theta$Omega1.inv
        } else {
          for (i in subj1) {
            Omega.inv.p2 <- Omega.inv.p2 + (theta$beta1[i, ]) %*% t(theta$beta1[i, ])
          }
          Omega.inv.mat <- solve(Lk + Omega.inv.p2)
          new.Omega1.inv <- MCMCpack::rwish(v = Omega.inv.df, S = Omega.inv.mat)
        }
        return(new.Omega1.inv)
      }
      
      proposal_function_Omega2_inv <- function() {
        ## Sample Omega^-1 from a Wishart distribution
        tau <- theta$tau
        n2 <- sum(tau == 1)
        subj2 <- which(tau == 1)
        new.Omega2.inv <- theta$Omega2.inv
        Omega.inv.df <- nu + n2 + 1
        Omega.inv.p2 <- matrix(0, nrow = K, ncol = K)
        if (length(subj2) == 0){
          new.Omega2.inv <- theta$Omega2.inv
        } else {
          for (i in subj2) {
            Omega.inv.p2 <- Omega.inv.p2 + (theta$beta2[i, ]) %*% t(theta$beta2[i, ])
          }
          Omega.inv.mat <- solve(Lk + Omega.inv.p2)
          new.Omega2.inv <- MCMCpack::rwish(v = Omega.inv.df, S = Omega.inv.mat)
        }
        return(new.Omega2.inv)
      }
      
      proposal_function_Zi1 <- function(update.all = FALSE) {
        new.Zi1 <- theta$Zi1
        tau <- theta$tau
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new.Zi1 <- theta$Zi1
        } else{
          if (update.all){
            for (i in subj1){
              Ai.temp <- diag(theta$D1[[id.all[i]]]) %*% theta$Sigma1.inv[[id.all[i]]]
              Ai <- diag(N.tps.all[i]) + Ai.temp %*% diag(theta$D1[[id.all[i]]])
              ai <- Ai.temp %*% (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta1[i, ])
              Zi.varcov <- solve(Ai)
              Zi.mean <- Zi.varcov %*% ai
              new.Zi1[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), 
                                                       mean = Zi.mean, sigma = Zi.varcov, 
                                                       lower = rep(0, N.tps.all[i]), 
                                                       upper = rep(Inf, N.tps.all[i]))$samp
            }
          } else{
            # Sample one value from subj1
            i <- sample(subj1, 1)
            Ai.temp <- diag(theta$D1[[id.all[i]]]) %*% theta$Sigma1.inv[[id.all[i]]]
            Ai <- diag(N.tps.all[i]) + Ai.temp %*% diag(theta$D1[[id.all[i]]])
            ai <- Ai.temp %*% (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta1[i, ])
            Zi.varcov <- solve(Ai)
            Zi.mean <- Zi.varcov %*% ai
            new.Zi1[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), 
                                                     mean = Zi.mean, sigma = Zi.varcov, 
                                                     lower = rep(0, N.tps.all[i]), 
                                                     upper = rep(Inf, N.tps.all[i]))$samp
          }
        }
        
        return(new.Zi1)
      }
      
      proposal_function_Zi2 <- function(update.all = FALSE) {
        new.Zi2 <- theta$Zi2
        tau <- theta$tau
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new.Zi2 <- theta$Zi2
        } else{
          if (update.all){
            for (i in subj2){
              Ai.temp <- diag(theta$D2[[id.all[i]]]) %*% theta$Sigma2.inv[[id.all[i]]]
              Ai <- diag(N.tps.all[i]) + Ai.temp %*% diag(theta$D2[[id.all[i]]])
              ai <- Ai.temp %*% (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta2[i, ])
              Zi.varcov <- solve(Ai)
              Zi.mean <- Zi.varcov %*% ai
              new.Zi2[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), 
                                                       mean = Zi.mean, sigma = Zi.varcov, 
                                                       lower = rep(0, N.tps.all[i]), 
                                                       upper = rep(Inf, N.tps.all[i]))$samp
            }
          } else{
            # Sample one value from subj2
            i <- sample(subj2, 1)
            Ai.temp <- diag(theta$D2[[id.all[i]]]) %*% theta$Sigma2.inv[[id.all[i]]]
            Ai <- diag(N.tps.all[i]) + Ai.temp %*% diag(theta$D2[[id.all[i]]])
            ai <- Ai.temp %*% (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta2[i, ])
            Zi.varcov <- solve(Ai)
            Zi.mean <- Zi.varcov %*% ai
            new.Zi2[[id.all[i]]] <- tmvnsim::tmvnsim(nsamp = 1, k = nrow(Zi.varcov), 
                                                     mean = Zi.mean, sigma = Zi.varcov, 
                                                     lower = rep(0, N.tps.all[i]), 
                                                     upper = rep(Inf, N.tps.all[i]))$samp
          }
        }
        
        return(new.Zi2)
      }
      
      proposal_function_D1 <- function() {
        new.D1 <- theta$D1
        tau <- theta$tau
        subj1 <- which(tau == 2)
        if (length(subj1) == 0){
          new.D1 <- theta$D1
        } else {
          for (i in subj1){
            B <- diag(x = 1/gamma, nrow = N.tps.all[i]) + diag(theta$Zi1[[id.all[i]]]) %*% 
              theta$Sigma1.inv[[id.all[i]]] %*% diag(theta$Zi1[[id.all[i]]])
            b <- diag(theta$Zi1[[id.all[i]]]) %*% theta$Sigma1.inv[[id.all[i]]] %*%
              (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta1[i, ])
            vecD.mean <- solve(B) %*% b
            vecD.varcov <- solve(B)
            new.D1[[id.all[i]]] <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
          }
        }
        return(new.D1)
      }
      
      proposal_function_D2 <- function() {
        new.D2 <- theta$D2
        tau <- theta$tau
        subj2 <- which(tau == 1)
        if (length(subj2) == 0){
          new.D2 <- theta$D2
        } else {
          for (i in subj2){
            B <- diag(x = 1/gamma, nrow = N.tps.all[i]) + diag(theta$Zi2[[id.all[i]]]) %*% 
              theta$Sigma2.inv[[id.all[i]]] %*% diag(theta$Zi2[[id.all[i]]])
            b <- diag(theta$Zi2[[id.all[i]]]) %*% theta$Sigma2.inv[[id.all[i]]] %*%
              (data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta2[i, ])
            vecD.mean <- solve(B) %*% b
            vecD.varcov <- solve(B)
            new.D2[[id.all[i]]] <- mnormt::rmnorm(n = 1, mean = vecD.mean, varcov = vecD.varcov)
          }
        }
        return(new.D2)
      }
      
      proposal_function_Sigma1_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        new.Sigma1.inv <- theta$Sigma1.inv
        tau <- theta$tau
        subj1 <- which(tau == 2)
        n1 <- sum(tau == 2)
        if (length(subj1) == 0){
          new.Sigma1.inv = theta$Sigma1.inv
        } else{
          for (i in subj1) {
            Sigma.inv.df <- twor.all[i] + n1
            Sigma.inv.mat <- solve(2 * kappa.list[[i]])
            temp <- data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta1[i, ] - diag(theta$D1[[id.all[i]]]) %*% theta$Zi1[[id.all[i]]]
            Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
            new.Sigma1.inv[[id.all[i]]] <- MCMCpack::rwish(v = Sigma.inv.df, S = solve(Sigma.inv.mat))
          }
        }
        return(new.Sigma1.inv)
      }
      
      proposal_function_Sigma2_inv <- function() {
        ## Sample Sigma^(-1) from a Wishart distribution
        new.Sigma2.inv <- theta$Sigma2.inv
        tau <- theta$tau
        subj2 <- which(tau == 1)
        n2 <- sum(tau == 1)
        if (length(subj2) == 0){
          new.Sigma2.inv = theta$Sigma2.inv
        } else{
          for (i in subj2) {
            Sigma.inv.df <- twor.all[i] + n2
            Sigma.inv.mat <- solve(2 * kappa.list[[i]])
            temp <- data$x.demean[[i]] - HUk.list[[i]] %*% theta$beta2[i, ] - diag(theta$D2[[id.all[i]]]) %*% theta$Zi2[[id.all[i]]]
            Sigma.inv.mat <- Sigma.inv.mat + temp %*% t(temp)
            new.Sigma2.inv[[id.all[i]]] <- MCMCpack::rwish(v = Sigma.inv.df, S = solve(Sigma.inv.mat))
          }
        }
        return(new.Sigma2.inv)
      }
      
      proposal_function_wi <- function() {
        new.wi <- theta$wi
        tau <- theta$tau
        ## sample wi from a Gamma distribution
        vi <- LaplacesDemon::rtrunc(N.subj, spec="gamma", a=2, b=Inf, shape = 1, rate = 0.1)
        for (i in 1:N.subj){
          if (tau[i] == 1){ #mvt
            wi.temp <- var(data$x.demean[[i]])
            wi.shape <- vi[i]/2 + N.subj/2
            wi.rate <- vi[i]/2 +  N.subj * wi.temp/ 2
            new.wi[i] <- rgamma(1, shape = wi.shape, rate = wi.rate)
          } else{ #mvn
            new.wi[i] <- 1
          }
        }
        return(new.wi)
      }
      
      proposal_function_pi1 <- function(){
        tau <- theta$tau
        n1 <- sum(tau == 2)
        n2 <- sum(tau == 1)
        new.pi1 <- rbeta(1, 1+n1, 1+n2)
        return(new.pi1)
      }
      
      proposal_function_tau <- function(){
        new.tau <- theta$tau
        tau <- theta$tau
        pi1 <- theta$pi1
        for (i in 1:N.subj){
          tau.A <- pi1 * LaplacesDemon::dmvn(x = data$x.demean[[i]],
                                             mu = t(HUk.list[[i]] %*% theta$beta1[i,] +
                                                      diag(theta$D1[[id.all[i]]]) %*% theta$Zi1[[id.all[i]]]),
                                             Sigma = solve(theta$Sigma1.inv[[id.all[i]]]))
          tau.B <- (1 - pi1) * LaplacesDemon::dmvn(x = data$x.demean[[i]],
                                                   mu = t(HUk.list[[i]] %*% theta$beta2[i,] +
                                                            diag(theta$D2[[id.all[i]]]) %*% theta$Zi2[[id.all[i]]]),
                                                   Sigma = solve(theta$Sigma2.inv[[id.all[i]]])/theta$wi[i])
          if (tau.A + tau.B == 0){
            new.tau[i] <- LaplacesDemon::rbern(1, 0) + 1
          } else{
            new.tau[i] <- LaplacesDemon::rbern(1, tau.A/(tau.A + tau.B)) + 1
          }
        }
        
        new.tau <- ReMap.mm(new.tau)
        
        return(new.tau)
      }
      
      # Map each parameter to its corresponding proposal function
      proposal_functions <- list(
        beta1 = proposal_function_beta1,
        Omega1.inv = proposal_function_Omega1_inv,
        Zi1 = proposal_function_Zi1,
        D1 =  proposal_function_D1,
        Sigma1.inv = proposal_function_Sigma1_inv,
        beta2 = proposal_function_beta2,
        Omega2.inv = proposal_function_Omega2_inv,
        Zi2 = proposal_function_Zi2,
        D2 =  proposal_function_D2,
        Sigma2.inv = proposal_function_Sigma2_inv,
        wi = proposal_function_wi,
        pi1 = proposal_function_pi1,
        tau = proposal_function_tau
      )
      
      # Update the chosen parameter using its respective proposal function
      if (chosen_param == "tau"){
        theta_new[["tau"]] <- proposal_functions[["tau"]]()
        theta_new[["beta1"]] <- proposal_functions[["beta1"]](update.all = TRUE)
        theta_new[["beta2"]] <- proposal_functions[["beta2"]](update.all = TRUE)
        theta_new[["Omega1.inv"]] <- proposal_functions[["Omega1.inv"]]()
        theta_new[["Omega2.inv"]] <- proposal_functions[["Omega2.inv"]]()
        theta_new[["Zi1"]] <- proposal_functions[["Zi1"]](update.all = TRUE)
        theta_new[["Zi2"]] <- proposal_functions[["Zi2"]](update.all = TRUE)
        theta_new[["D1"]] <- proposal_functions[["D1"]]()
        theta_new[["D2"]] <- proposal_functions[["D2"]]()
        theta_new[["Sigma1.inv"]] <- proposal_functions[["Sigma1.inv"]]()
        theta_new[["Sigma2.inv"]] <- proposal_functions[["Sigma2.inv"]]()
        theta_new[["wi"]] <- proposal_functions[["wi"]]()
        theta_new[["pi1"]] <- proposal_functions[["pi1"]]()
      } else{
        theta_new[[chosen_param]] <- proposal_functions[[chosen_param]]()
      }
      
      # compute numerator and denominator of the MH ratio
      num <- (likelihood(theta_new)) * alpha[[r]] +  prior(theta_new)
      den <- (likelihood(theta)) * alpha[[r]] + prior(theta)
      
      # compute ratio
      ratio <- min(1, exp(num - den))
      
      # accept/reject step
      if(runif(1) < ratio){
        theta = theta_new
      } else{
        theta = theta
      }
      
      return(theta)
    }
    
    particles <- parLapply(cl, previous_particles, Gibbs_HM_Move)
    
    # compute the ESS
    log_incremental_w <- alphaDiff * u
    logw <- log_incremental_w + logw  # log unnormalized weights
    logmax <- max(logw)
    logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  
    W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
    logW <- log(W[[r]])                                # log normalized weights
    cat("logZ:",logZ[[r]],"\n")
    ESS[[r]] <- rESS(logW)
    
    # resample if ESS below threshold
    if( ESS[[r]]<tuning_param$eps )
    {
      cat("Resample: ESS=", ESS[[r]], '\n')
      ancestors <- systematic_resample( W[[r]] )
      particles <-  particles[ancestors]
      W[[r]] <- rep(1/N.particles,N.particles)
      logW <- log(W[[r]])
      w <- rep(1,N.particles)
      logw <- rep(0,N.particles)
    }
    previous_particles <- particles
    
  }
  stopCluster(cl)
  
  return(list(particles = particles,
              alpha = alpha,
              ESS = ESS,
              logZ = logZ,
              W = W,
              H = H.list, Uk = Uk.list, Lk = Lk.list))
}
