sim_splines <- function(j,n,tau,cp,lam0,beta0,
                        xies_start,gammas_start,
                        tol=1e-5, maxit=1000){
  tryCatch( 
    { 
      source("functions/dgps_splines.R", local=TRUE)
      source("functions/em_splines.R", local=TRUE)
      set.seed(j)
      gc()
      # Generate data
      ones <- rep(1,n)
      age <- rpois(n,45)+20
      period1 <- runif(floor(n/3))*15*365
      period2 <- (15+runif(n-floor(n/3))*15)*365
      period <- c(period1,period2)
      sex <- as.numeric(runif(n)<0.48)
      socio_economic <- matrix(c(.3,.6,.1),n,3, byrow=TRUE) # low, mid, high
      socio_economic <- apply(socio_economic, 1, function(p) 
      sample.int(3L, size = 1L, prob = p ))
      socio_economic_dummy <- t(vapply(socio_economic, function(k){
        vec <- integer(3L) 
        vec[k] <- 1L 
        vec
      }, integer(3L)))
      tumor_size <- rchisq(n, df=2) # colon cancer average tumor size
      log_tms <- log(tumor_size)
      log_tms_s <- as.numeric(scale(log_tms))
      scaled_age <- scale(age)
      scaled_period <- scale(period)
      X <- as.data.frame(cbind(ones,scaled_age,scaled_period,sex, socio_economic_dummy, log_tms_s))
      colnames(X) <- c("ones","age", "period", "sex", "socio_low", "socio_mid", "socio_high", "log_tms")
      pies_true <- c(
        -0.405,  # intercept ~ logit(0.4)
        -0.30,   # age
        0.05,   # sex
        0.20,   # age^2
        0.10,   # period^2
        0.15,   # |period|*age
        -0.10,   # socio_low
        0.01,   # socio_mid
        0.10,   # socio_high
        -0.05    # log_tms*age
      )
      gammas_true <- c(-1.8,
                       0.001,
                       -0.2,
                       -0.001)
      #
      xies_true <- c(-1.8,
                     0.001,
                     -0.1,
                     -0.001)
      #
      X <- as.matrix(X)
      # Generate Data
      probs <- plogis(pies_true[1]*(X[,1]) +
                        pies_true[2]*(X[,2]) +
                        pies_true[3]*(X[,4]) +
                        pies_true[4]*(X[,2]^2) +
                        pies_true[5]*(X[,3]^2) +
                        pies_true[6]*(abs(X[,3])*X[,2]) +
                        pies_true[7]*(X[,5]) +
                        pies_true[8]*(X[,6]) +
                        pies_true[9]*(X[,7]) +
                        pies_true[10]*(log(tumor_size)*X[,2]))
      lamd <- exp(X[,1:4]%*%gammas_true)
      betad <- exp(X[,1:4]%*%xies_true)
      #
      sim <-dgps_splines(n, X, tau, probs, lamd, betad, lam0, beta0, cp)
      gc()
      res<-em_splines(t_obs = sim$t_tot,
                      deltas = sim$deltas,
                      X=X,
                      lam0 = lam0,
                      beta0 = beta0,
                      xies_start = xies_start,
                      gammas_start = gammas_start,
                      tol = tol,
                      maxit=maxit,
                      j=j)
      gc()
      cat("\n",j, "Parameter estimated", sep = " ")
      cat(j,
          mean(abs(res$prob0-probs)),
          mean((res$prob0-probs)^2),
          gammas_true,
          res$gammas,
          res$gammas-gammas_true,
          xies_true,
          res$xies,
          res$xies-xies_true,
          abs(tail(res$lliks, n=1)[1]-tail(res$lliks, n=2)[1]),
          ll_optim = tail(res$lliks, n=1),
          res$iter, 
          file="log/log.txt", sep = " ; ", append = TRUE)
      #
      write("\n",file="log/log.txt",append=TRUE)
      gc()
      #
      return(list(prob0_MAE = mean(abs(res$prob0-probs)),
                  prob0_MSE = mean((res$prob0-probs)^2),
                  gammas_true=gammas_true,
                  gammas = res$gammas,
                  gammas_diff = res$gammas-gammas_true,
                  xies_true=xies_true,
                  xies = res$xies,
                  xies_diff = res$xies-xies_true,
                  diff_last = (tail(res$lliks, n=1)[1]-tail(res$lliks, n=2)[1]),
                  ll_optim = tail(res$lliks, n=1),
                  iter=res$iter))
    },
    error = function(w){
      cat("Error at iteration:",j, 
          file="log/log.txt", sep = " ", append = TRUE)
      write("\n",file="log/log.txt",append=TRUE)
      cat("Error:",warning(w), 
          file="log/log.txt", sep = " ", append = TRUE)
      write("\n",file="log/log.txt",append=TRUE)
    }
  )
  #
}
