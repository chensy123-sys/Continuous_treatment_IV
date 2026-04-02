rm(list = ls());source('cross_fitting2.R')
n <- 5000; Looptime <- 200; taeget_point = 0.5; interval = c(0.25,0.75); test =NULL
theta <- function(A,U,L){return(A *1 + U- L*0.5)}
set.seed(2025);delta <- delta_adap_func(sim_data(1e4, theta_func = theta), target_point = taeget_point)

set.seed(2025);LLKR_fit_list <- list()
for (loop in 1:Looptime) {
  mydat <- sim_data(n,theta_func = theta)
  fit_CF <- cross_fitting(data = mydat,IV_func = delta,spline_k=30,K=5,constriant = 'L',print_state = T)
  # test <- URWF_test(data = mydat,func = delta,region = seq(-0.9,0.9,length.out=50),h = NULL,constriant = 'L')
  cl <- clean_extreme_scores(fit_CF, test = test, alpha = 0.01, interval = interval)
  LLKR_fit <- LLKR(cl); LLKR_fit_list[[loop]] <- LLKR_fit
  cat('Loop:',loop,'\n')
}
save(LLKR_fit_list, file = "LLKR_fit_2.RData")
