rm(list = ls()) # Clear all objects in the workspace
library(np)    # Nonparametric kernel methods
library(dplyr) # Data manipulation
source('cross_fitting2.R')
# Create temporary environments to load multiple RData files separately
env1 <- new.env()
env2 <- new.env()

load("data/LLKR_fit_1.RData", envir = env1)
load("data/LLKR_fit_2.RData", envir = env2)

# Combine LLKR_fit_list from both environments
LLKR_fit_list <- c(env1$LLKR_fit_list, env2$LLKR_fit_list)
Looptime <- length(LLKR_fit_list)  # Number of models
A_test <- c(0.4, 0.5, 0.6)      # Fixed treatment values for prediction

set.seed(2025)

# Initialize matrices to store predicted outcomes
Y_AIPW <- Y_IPW <- Y_OR <- matrix(0, Looptime, length(A_test))
Y_AIPW_nuc <- Y_IPW_nuc <- Y_OR_nuc <- matrix(0, Looptime, length(A_test))

# Loop through all LLKR fits and generate predictions for A_test
for(i in 1:Looptime){
  Y_AIPW[i,] <- predict(LLKR_fit_list[[i]]$LLKR_AIPW, newdata = data.frame(treatment=A_test))
  Y_IPW[i,] <- predict(LLKR_fit_list[[i]]$LLKR_IPW, newdata = data.frame(treatment=A_test))
  Y_OR[i,] <- predict(LLKR_fit_list[[i]]$LLKR_OR, newdata = data.frame(treatment=A_test))
  
  Y_AIPW_nuc[i,] <- predict(LLKR_fit_list[[i]]$LLKR_AIPW_nuc, newdata = data.frame(treatment=A_test))
  Y_IPW_nuc[i,] <- predict(LLKR_fit_list[[i]]$LLKR_IPW_nuc, newdata = data.frame(treatment=A_test))
  Y_OR_nuc[i,] <- predict(LLKR_fit_list[[i]]$LLKR_OR_nuc, newdata = data.frame(treatment=A_test))
  
  if(i %% 100 == 0) cat(i, ' ')
}

# Generate a larger set of treatment values using simulation
set.seed(2025)
A_test_2 <- sim_data(1e4, theta_func=theta) %>% select(A) %>% filter(A >= 0.25 & A <= 0.75)
A_test_2 <- A_test_2$A[1:1000]

# Initialize matrices for new predictions
Y_AIPW_2 <- Y_IPW_2 <- Y_OR_2 <- matrix(0, Looptime, length(A_test_2))
Y_AIPW_nuc_2 <- Y_IPW_nuc_2 <- Y_OR_nuc_2 <- matrix(0, Looptime, length(A_test_2))

# Loop through all LLKR fits for the second test set
for(i in 1:Looptime){
  Y_AIPW_2[i,] <- predict(LLKR_fit_list[[i]]$LLKR_AIPW, newdata = data.frame(treatment=A_test_2))
  Y_IPW_2[i,] <- predict(LLKR_fit_list[[i]]$LLKR_IPW, newdata = data.frame(treatment=A_test_2))
  Y_OR_2[i,] <- predict(LLKR_fit_list[[i]]$LLKR_OR, newdata = data.frame(treatment=A_test_2))
  
  Y_AIPW_nuc_2[i,] <- predict(LLKR_fit_list[[i]]$LLKR_AIPW_nuc, newdata = data.frame(treatment=A_test_2))
  Y_IPW_nuc_2[i,] <- predict(LLKR_fit_list[[i]]$LLKR_IPW_nuc, newdata = data.frame(treatment=A_test_2))
  Y_OR_nuc_2[i,] <- predict(LLKR_fit_list[[i]]$LLKR_OR_nuc, newdata = data.frame(treatment=A_test_2))
  
  if(i %% 100 == 0) cat(i, ' ')
}

my_tab <- data.frame(matrix(0,6,8), row.names = c('IV_AIPW','IV_IPW','IV_OR','NUC_AIPW','NUC_IPW','NUC_OR'))
colnames(my_tab) <- c('bias_0.4','bias_0.5','bias_0.6','bias_N', 'rmse_0.4','rmse_0.5','rmse_0.6','rmse_N')
temp <- matrix(rep(A_test, Looptime), Looptime, length(A_test), byrow = TRUE)
temp_2 <- matrix(rep(A_test_2, Looptime), Looptime, length(A_test_2), byrow = TRUE)
c(round(abs(apply(Y_AIPW, 2, mean) - A_test), 4),
  round(mean(abs(apply(Y_AIPW_2, 2, mean) - A_test_2)), 4),
  round(sqrt(apply((Y_AIPW-temp)^2,2,mean)) ,4),
  round(sqrt(mean(apply((Y_AIPW_2 - temp_2)^2, 2, mean))), 4)) -> my_tab[1,]

c(round(abs(apply(Y_IPW, 2, mean) - A_test), 4),
  round(mean(abs(apply(Y_IPW_2, 2, mean) - A_test_2)), 4),
  round(sqrt(apply((Y_IPW-temp)^2,2,mean)) ,4),
  round(sqrt(mean(apply((Y_IPW_2 - temp_2)^2, 2, mean))), 4)) -> my_tab[2,]

c(round(abs(apply(Y_OR, 2, mean) - A_test), 4),
  round(mean(abs(apply(Y_OR_2, 2, mean) - A_test_2)), 4),
  round(sqrt(apply((Y_OR-temp)^2,2,mean)) ,4),
  round(sqrt(mean(apply((Y_OR_2 - temp_2)^2, 2, mean))), 4)) -> my_tab[3,]



c(round(abs(apply(Y_AIPW_nuc, 2, mean) - A_test), 4),
  round(mean(abs(apply(Y_AIPW_nuc_2, 2, mean) - A_test_2)), 4),
  round(sqrt(apply((Y_AIPW_nuc-temp)^2,2,mean)) ,4),
  round(sqrt(mean(apply((Y_AIPW_nuc_2 - temp_2)^2, 2, mean))), 4)) -> my_tab[4,]

c(round(abs(apply(Y_IPW_nuc, 2, mean) - A_test), 4),
  round(mean(abs(apply(Y_IPW_nuc_2, 2, mean) - A_test_2)), 4),
  round(sqrt(apply((Y_IPW_nuc-temp)^2,2,mean)) ,4),
  round(sqrt(mean(apply((Y_IPW_nuc_2 - temp_2)^2, 2, mean))), 4)) -> my_tab[5,]

c(round(abs(apply(Y_OR_nuc, 2, mean) - A_test), 4),
  round(mean(abs(apply(Y_OR_nuc_2, 2, mean) - A_test_2)), 4),
  round(sqrt(apply((Y_OR_nuc-temp)^2,2,mean)) ,4),
  round(sqrt(mean(apply((Y_OR_nuc_2 - temp_2)^2, 2, mean))), 4)) -> my_tab[6,]


write.csv(my_tab, file = 'data/TableC.1(n=5000).csv')
