dgps_nn<-function(n, X, tau, probs_vec, lamd_vec, betad, lam0, beta0, cp, maxtry=1000){
  
  # n: number of observations
  # X: matrix of covariates' values for the n observations and an initial column of ones
  # tau: maximum follow-up time in days
  # probs_vec: vector of probability of being cured for each individual
  # lamd_vec: vector of excess hazard of each individual
  # betad: Shape excess surv
  # lam0: Scale baseline/expected surv
  # beta0: Shape baseline/expected surv
  # cp: censoring proportion
  
  # z=1 cured, z=0 uncured (convention baseline surv S0)
  z <- as.numeric(runif(n)<probs_vec)
  n1 <- sum(z)
  n0 <- n-n1
  # Generate time to event for cured proportion 
  u <- runif(n1) # u ~ U[0,1] iff 1-u ~ U[0,1]
  t_cure <- ((-log(u))^(1/beta0))/lam0
  # Generate time to event for uncured proportion 
  u <- runif(n0)
  finvuc <- function(x, u, lamd, betad){
    u-exp(-(lam0*x)^(beta0)-(lamd*x)^(betad))
  }
  t_uncure <- rep(NA, n0)
  lamd_sel <- lamd_vec[as.logical(z)] # Select unit-specific hazard
  betad_sel <- betad[as.logical(z)] # Select unit-specific hazard
  env_dgps<-new.env()
  assign("go", FALSE, env=env_dgps)
  assign("ntry", 1, env=env_dgps)
  assign("x", 0, env=env_dgps)
  for (i in 1:n0){
    assign("unf", u[i], env=env_dgps)
    while ((!get("go", env=env_dgps))&(get("ntry", env=env_dgps)<maxtry)) {
      tryCatch(
        {
          assign("x",uniroot(finvuc, interval = c(0,tau),
                             extendInt = "upX", tol=1e-12,
                             u = get("unf", env=env_dgps),
                             lamd = lamd_sel[i], betad=betad_sel[i])$root,
                 env=env_dgps)
          assign("go", TRUE, env=env_dgps)
          
        },
        error = function(e){
          assign("unf", runif(1), env=env_dgps)
          assign("ntry", get("ntry", env=env_dgps)+1, env=env_dgps)
        }
      )
    }
    t_uncure[i] <- get("x", env=env_dgps)
    assign("go", FALSE, env=env_dgps)
    assign("ntry", 1, env=env_dgps)
  }
  t_uncure[t_uncure<1e-32] <- 1e-32
  t_tot <- rep(NA,n)
  t_tot[z==1]<-t_cure
  t_tot[z==0]<-t_uncure
  # Censoring
  deltas <- rep(1,n)
  deltas[which(t_tot>tau)] <- 0
  t_tot[deltas==0] <- tau
  t_tot[t_tot<1e-32]<-1e-32
  return(list(t_tot=t_tot,deltas=deltas,z=z))
} 