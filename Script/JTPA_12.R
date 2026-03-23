library(dplyr);library(tibble);library(tidyr);library(readr)
set.seed(2025)
rm(list = ls()); source('cross_fitting2.R')
mydat <- read.csv("jpta_han.csv") %>% select(-X) %>% 
  mutate(
    A = pmax(pmin(A,15),8),
    # Y = as.numeric(Y>median(Y)),
  )
target_point = 12; constriant = "L"; spline_k =4

delta <- delta_adap_func(mydat, target_point = target_point,constriant = NULL,h = 1)
Looptime <- 400
fit_CF_list = list()
for (loop in 1:Looptime) {
  mydat_boot <- mydat[sample(nrow(mydat), nrow(mydat), replace = TRUE),]
  
  fit_CF_list[[loop]] <- cross_fitting(
    data = mydat_boot,
    IV_func = delta,
    spline_k = spline_k, print_state = T,
    K=5,constriant = constriant
  )  
  print(loop)
}

save.image(file = paste0('data/JTPA_',target_point,'.RDATA'))