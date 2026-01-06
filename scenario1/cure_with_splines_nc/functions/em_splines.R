em_splines <- function(t_obs, deltas, X, lam0, beta0, xies_start, gammas_start, tol = 1e-6,maxit=1000, j=1){
  #
  env<-new.env()
  n <- length(t_obs)
  pp1 <- dim(X)[2] # p + 1  ; covariates + column of ones
  #
  h0 <- beta0*lam0*(lam0*t_obs)^(beta0-1) # baseline hazard
  S0 <- exp(-(lam0*t_obs)^(beta0)) # baseline survival
  #
  log_lik <- function(deltas, prob0, S0, h0, Sd, hd){
    fab <- prob0*S0*(h0^deltas) + (1-prob0)*S0*Sd*((h0 + hd)^(deltas))
    fab[fab<1e-32] <- 1e-32
    sum(log(fab))
  }
  #
  complete_loglik <- function(params, lb){
    ls <- exp(X%*%lb[1:params$pp1])
    lbetad <- exp(X%*%lb[(params$pp1+1):(2*params$pp1)]) 
    fab <- exp(-((ls*params$t_obs)^lbetad))*
      ((params$h0+lbetad*ls*((ls*params$t_obs)^(lbetad-1)))^params$deltas)
    fab[fab<1e-32] <- 1e-32
    sum((1-params$z)*log(fab))
  }
  #
  z <- rep(NA, n)
  assign("gammas", gammas_start, env=env)
  assign("xies", xies_start, env=env)
  prob0 <- rep(0.5, n)
  lamd <- exp(X%*%get("gammas", env=env))
  betad <- exp(X%*%get("xies", env=env))
  #
  hd <- betad*lamd*((lamd*t_obs)^(betad-1)) # Excess hazard
  Sd <- exp(-(lamd*t_obs)^betad) # Excess survival
  #
  y1 <- rep(1,n) # cured
  y0 <- rep(0,n) # uncured
  y <- c(y1,y0)
  #
  # standardize features and construct splines bases on age and period
  X_gam <- as.data.frame(rbind(X,X))
  #
  lliks<-c()
  log_lik_old <- log_lik(deltas, prob0, S0, h0, Sd, hd)
  iter <- 0
  dif <- 1
  cat("iter: ", 0," dif: ",dif," llk: ", log_lik_old,"\n")
  #
  while ((dif > tol) & (iter<maxit)) { 
    # Expectation
    fab <- prob0*S0*(h0^deltas)
    z <- fab/((1-prob0)*S0*Sd*((h0+hd)^deltas)+fab)
    z[z<1e-32]<-1e-32
    #
    # Maximization
    # Update pies
    w <- c(z,1-z)
    tryCatch(
      {
        assign("gam_fit", gam(y ~ sex + s(age)+s(period)+ti(age,period), 
                              family = binomial(), 
                              data=X_gam, weights = w) , env=env)
      }, error=function(e){
        assign("gam_fit", gam(y ~ sex + s(age)+s(period)+ti(age,period), 
                              family = binomial(), 
                              data=X_gam, weights = w, method = "REML") , env=env)
      }
    )
    gam_fit <- get("gam_fit", env=env)
    prob0 <- gam_fit$fitted.values[1:n]
    #
    assign("go", FALSE, env=env)
    assign("ntry", 0, env=env)
    while ((!get("go", env=env))&(get("ntry", env=env)<maxit)) {
      tryCatch(
        {
          assign("res",optim(par = c(get("gammas", env=env),get("xies", env=env)), 
                             fn = complete_loglik,
                             method = "Nelder-Mead",
                             params = list(t_obs=t_obs, z=z, deltas=deltas,
                                           S0=S0, h0=h0, pp1 = pp1),
                             control=list(fnscale=-1)),env=env)
          assign("go", TRUE, env=env)
        },
        error = function(e){
          if(get("ntry", env=env)==1){
            assign("gammas",gammas_start + c(runif(1, min =  0.0001,max = 0.001),
                                             runif(1, min =  0.0001,max = 0.001), 
                                             runif(1, min =  0.0001,max = 0.001),
                                             runif(1, min =  0.0001,max = 0.001)),
                   env = env)
            assign("xies",xies_start + c(runif(1, min =  0.0001,max = 0.001),
                                         runif(1, min =  0.0001,max = 0.001), 
                                         runif(1, min =  0.0001,max = 0.001),
                                         runif(1, min =  0.0001,max = 0.001)),
                   env = env)
          }
          else{
            assign("gammas", get("gammas", env=env) + c(runif(1, min =  0.0001,max = 0.001),
                                                        runif(1, min =  0.0001,max = 0.001), 
                                                        runif(1, min =  0.0001,max = 0.001),
                                                        runif(1, min =  0.0001,max = 0.001)), 
                   env=env)
            assign("xies", get("xies", env=env) + c(runif(1, min =  0.0001,max = 0.001),
                                                    runif(1, min =  0.0001,max = 0.001), 
                                                    runif(1, min =  0.0001,max = 0.001),
                                                    runif(1, min =  0.0001,max = 0.001)), 
                   env=env) 
          }
          assign("ntry", get("ntry", env=env)+1, env=env)
        }
      )
    }
    res <- get("res", env=env)
    assign("gammas", res$par[1:pp1], env=env)
    assign("xies", res$par[(pp1+1):(2*pp1)], env=env)
    #
    #
    # ...........................
    #
    # log-likelihood calculation:
    lamd <- exp(X%*%get("gammas", env=env))
    betad <- exp(X%*%get("xies", env=env))
    Sd <- exp(-(lamd*t_obs)^(betad)) # excess survival
    hd <- betad*lamd*(lamd*t_obs)^(betad-1) # excess hazard
    #
    log_lik_cur <- log_lik(deltas, prob0, S0, h0, Sd, hd)
    #
    # prepare for next iteration
    iter <- iter + 1 
    lliks <- c(lliks,log_lik_cur)
    dif <-  abs(log_lik_cur - log_lik_old)
    cat(j,"-", "iter: ", iter," dif: ", dif," llk: ",log_lik_cur,"\n")
    log_lik_old <- log_lik_cur
    gc()
  }
  return(list(lliks=lliks,prob0=prob0,gammas=get("gammas", env=env),xies=get("xies", env=env), iter=iter))
} 