require(parallel)
source("functions/sim_glm.R")
# Common parameters setting ----
tau <- 30*365 # 30 years follow-up: the follow up time in the data is in days!
lam0 <- 0.00009 # Scale baseline/expected surv
beta0 <- 30 # Shape baseline/expected surv
cp <- 0.08 # censoring proportion:
nsim <- 250  # number of simulation
# Simulation study ---- 
i<-1
for (n in c(10000,1000,500,100)){  
  if (n == 10000){
    cons_check <- as.data.frame(matrix(NA,nrow = 4, ncol = 18))
    row.names(cons_check)<-c("10000","1000","500","100")
    colnames(cons_check) <- c("mean_MAE_prob0",
                              "bias_gamma1",
                              "bias_gamma2",
                              "bias_gamma3",
                              "bias_gamma4",
                              "bias_xies1",
                              "bias_xies2",
                              "bias_xies3",
                              "bias_xies4",
                              "mean_MSE_prob0",
                              "MSE_gamma1",
                              "MSE_gamma2",
                              "MSE_gamma3",
                              "MSE_gamma4",
                              "MSE_xies1",
                              "MSE_xies2",
                              "MSE_xies3",
                              "MSE_xies4")
  }
  cl <- makeCluster(16, outfile="log/log_workers.txt")
  res <- parLapply(cl,1:nsim,sim_glm,
                   n,tau,cp,lam0,beta0,
                   gammas_start = c(0.0001,0.0001,0.0001,0.0001),
                   xies_start = c(0.0001,0.0001,0.0001,0.0001),
                   maxit=100)
  stopCluster(cl)
  gc()
  assign(paste("Res",n,sep = ""), matrix(unlist(res), nrow = nsim, ncol = 29,byrow=TRUE),env =.GlobalEnv)
  fab <- get(paste("Res",n,sep = ""))
  colnames(fab)<-c("prob0_MAE","prob0_MSE",
                   "gamma1_true","gamma2_true","gamma3_true","gamma4_true",
                   "gamma1","gamma2","gamma3","gamma4",
                   "gamma1_diff","gamma2_diff","gamma3_diff","gamma4_diff",
                   "xies1_true","xies2_true","xies3_true","xies4_true",
                   "xies1","xies2","xies3","xies4",
                   "xies1_diff","xies2_diff","xies3_diff","xies4_diff",
                   "diff_last","ll_optim","iter")
  assign(paste("Res",n,sep = ""), fab, env =.GlobalEnv)
  ## Check consistency of estimators ----
  cons_check[i,] <- c(apply(fab,2,mean)[1],  # mean_MAE_prob0
                      apply(fab,2,mean)[11], # Bias gamma1
                      apply(fab,2,mean)[12], # Bias gamma2
                      apply(fab,2,mean)[13], # Bias gamma3
                      apply(fab,2,mean)[14], # Bias gamma4
                      apply(fab,2,mean)[23], # Bias xi1
                      apply(fab,2,mean)[24], # Bias xi2
                      apply(fab,2,mean)[25], # Bias xi3
                      apply(fab,2,mean)[26], # Bias xi4
                      apply(fab,2,mean)[2], # mean_MSE_prob0
                      apply(fab^2,2,mean)[11], # MSE gamma1
                      apply(fab^2,2,mean)[12], # MSE gamma2
                      apply(fab^2,2,mean)[13], # MSE gamma3
                      apply(fab^2,2,mean)[14], # MSE gamma4
                      apply(fab^2,2,mean)[23], # MSE xi1
                      apply(fab^2,2,mean)[24], # MSE xi2
                      apply(fab^2,2,mean)[25], # MSE xi3
                      apply(fab^2,2,mean)[26]) # MSE xi4
  i<-i+1
  cat("End iters for n: ", n,"\n", file = "log/log.txt", append = TRUE)
}
#
save.image(file="data/data_glm.RData") 

