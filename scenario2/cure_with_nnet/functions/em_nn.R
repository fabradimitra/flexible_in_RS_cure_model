em_nn <- function(t_obs, deltas, X, lam0, beta0, xies_start, gammas_start, size, tol = 1e-6,maxit=1000, j){
  #
  neg_dif <- 0
  env<-new.env()
  n <- length(t_obs)
  pp1 <- 4 # only age period and sex
  pp1_nnet <- dim(X)[2] # p + 1  ; covariates + column of ones
  X_nn <- X 
  X <- X_nn[,1:pp1] # only
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
    ls <- exp(params$X%*%lb[1:params$pp1])
    lbetad <- exp(params$X%*%lb[(params$pp1+1):(2*params$pp1)]) 
    fab <- exp(-((ls*params$t_obs)^lbetad))*
      ((params$h0+lbetad*ls*((ls*params$t_obs)^(lbetad-1)))^params$deltas)
    fab[fab<1e-32] <- 1e-32
    sum((1-params$z)*log(fab))
  }
  #
  nn_cure_loss <- function(z,p){
    -sum(z*log(p) + (1-z)*log(1-p))
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
  X_nnet <- as.matrix(rbind(X_nn,X_nn))
  y0 <- rep(0,n)
  y1 <- rep(1,n)
  y <- c(y1,y0)
  #
  lliks<-c()
  log_lik_old <- log_lik(deltas, prob0, S0, h0, Sd, hd)
  iter <- 0
  dif <- 1
  cat(j,"-", "iter: ", 0," dif: ",dif," llk: ", log_lik_old,"\n")
  #
  while ((dif > tol) & (iter<maxit)) {
    # Expectation
    fab <- prob0*S0*(h0^deltas)
    z <- fab/((1-prob0)*S0*Sd*((h0+hd)^deltas)+fab)
    z[z<1e-32]<-1e-32
    #
    w <- c(z,1-z)
    # Maximization
    # Train neural network
    if(iter>0){
      loss_old <- nn_cure_loss(z,prob0)
      fab_weights <- nn_cure$wts
      nn_cure <- nnet(x=X_nnet[,2:pp1_nnet], y=y, Wts=fab_weights,
                      weights = w,
                      trace=FALSE, Hess=FALSE,
                      size=size, maxit = 100, entropy = TRUE) # link logit
    }else{
      loss_old <- Inf
      set.seed(123)
      nn_cure <- nnet(x=X_nnet[,2:pp1_nnet], y=y, weights = w,
                      trace=FALSE, Hess=FALSE,
                      size=size, maxit = 1, entropy = TRUE) # link logit
    }
    #
    prob0 <- predict(nn_cure,X_nn[,2:pp1_nnet])
    #
    # Update gammas and beta
    assign("go", FALSE, env=env)
    assign("ntry", 0, env=env)
    while ((!get("go", env=env))&(get("ntry", env=env)<maxit)) {
      tryCatch(
        {
          assign("res",optim(par = c(get("gammas", env=env),get("xies", env=env)), 
                             fn = complete_loglik,
                             method = "Nelder-Mead",
                             params = list(t_obs=t_obs, z=z, deltas=deltas,
                                           S0=S0, h0=h0, pp1 = pp1, X=X),
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
    dif <-  log_lik_cur - log_lik_old
    if(dif<0){
      neg_dif <- neg_dif + 1
    }
    dif <- abs(dif)
    cat(j,"-", "iter: ", iter," dif: ", dif," llk: ",log_lik_cur,"\n")
    log_lik_old <- log_lik_cur
  }
  return(list(lliks=lliks,prob0=prob0,gammas=get("gammas", env=env),xies=get("xies", env=env), iter=iter, neg_dif = neg_dif))
} 