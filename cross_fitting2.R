library(dplyr);library(np);library(mgcv)
theta <- function(A,U,L){return(A *1 + U- L*0.5)}
sim_data <- function(n, theta_func, A_value = NULL){
  L=(runif(n)-0.5)
  U = (runif(n)-0.5)*3
  Z = -0.5*L + (runif(n)-0.5)*3
  sel = rbinom(n,1,0.7)
  if(is.null(A_value)){
    a <- plogis(Z*2+rnorm(n)) 
    b <- plogis(-U*2+rnorm(n))
    c <- 40
    A = 2 * ifelse(
      sel,
      a,#rbeta(n, c*a, c*(1-a)),
      b#rbeta(n, c*b, c*(1-b))
    ) - 1
  }else{
    A = rep(A_value,n)
  }
  res <-theta_func(A,U,L)
  Y = res
  return(data.frame(
    L=L,Z=Z,A=A,Y=Y,sel=sel
  ))  
}

sim_data_set_A <- function(theta_func,A_value,n=1e5){
  res <- rep(0,length(A_value))
  for(i in 1:length(A_value)){
    mydat <- sim_data(n,theta_func = theta_func, A_value = A_value[i])
    res[i] <- mean(mydat$Y)
  }
  return(res)
}


#' ------------------------------------------------------------------------------
#' Adaptive regular weighting function using a binary neighborhood indicator
#'
#' Description:
#'   Estimates the ratio
#'       P(|A - a0| <= 1 | Z, L) / P(|A - a0| <= 1 | L)
#'   using two GAM models. This ratio is later used as an adaptive
#'   instrument weighting function.
#'
#' Arguments:
#'   data         : data.frame containing variables A, Z, L
#'   target_point : target value a0 of treatment A
#'   h           : bandwidth 
#'
#' Returns:
#'   A function that takes newdata and returns the estimated ratio
#' ------------------------------------------------------------------------------
delta_adap_func <- function(data, target_point, h = NULL, constriant = NULL) {
  
  A <- data$A; Z <- data$Z
  
  # Indicator for a neighborhood around the target treatment value
  if(is.null(h)){
    h <- 1.06*sd(A) * length(A)^(-1/5)  # Silverman rule of thumb
  }
  # Epanechnikov kernel weights
  u <- (A - target_point) / h
  K_weights <- 0.75 * (1 - u^2) * (abs(u) <= 1)
  
  
  if(!is.null(constriant)){
    # L is not empty
    L <- data[[constriant]]
    m <- length(unique(L))
    if(m < 5){
      # L is discrete
      mean_model <- gam(K_weights ~ s(Z, by = factor(L)) + factor(L), method = "REML")
      mean_model2 <- gam(K_weights ~ factor(L), method = "REML")    
    }else{
      # L is continuous
      mean_model <- gam(K_weights ~ s(Z, L, k = 40), method = "REML")
      mean_model2 <- gam(K_weights ~ s(L, k = 40), method = "REML")
    }
  }else{
    # L is empty
    mean_model <- gam(K_weights ~ s(Z), method = "REML")
    mean_model2 <- gam(K_weights ~ 1, method = "REML")
  }
  
  
  # Return prediction function
  function(newdata) {
    pred1 <- predict(mean_model, newdata = newdata, type = "response")
    pred2 <- predict(mean_model2, newdata = newdata, type = "response")
    pred1 / pred2
  }
}


#' -------------------------------------------------------------------------
#' Function: cross_fitting
#'
#' Description:
#'   Perform K-fold cross-fitting to estimate several causal scores:
#'     - OR  (Outcome Regression)
#'     - IPW (Inverse Probability Weighting style score)
#'     - AIPW (Augmented IPW / Doubly Robust score)
#'
#'   The function also computes "nuisance-unadjusted" versions of these scores.
#'
#' Instrument construction:
#'   - If IV_func is a function: it is treated as a pre-specified RWF (random
#'     weighting function) and directly applied to the data.
#'   - If IV_func is a numeric value: an adaptive RWF is trained using
#'     delta_adap_func at the specified target point.
#'
#' Cross-fitting structure:
#'   - The dataset is split into K folds.
#'   - For each fold:
#'       * nuisance functions are estimated on the training data
#'       * scores are evaluated on the validation data
#'
#' Inputs:
#'   data        : data.frame containing variables A, Y, Z, L (required columns)
#'   K           : number of cross-fitting folds (default = 5)
#'   proportion  : proportion of data used for adaptive RWF training
#'   IV_func     : either
#'                   (1) a function mapping data -> instrument value
#'                   (2) a numeric target point used to train an adaptive RWF
#'   print_state : if TRUE, print fold progress
#'   spline_k    : number of spline basis functions used in GAM models
#'
#' Outputs:
#'   A list containing:
#'     $score : data.frame with estimated scores evaluated at each A
#'     $K     : number of folds used
#' -------------------------------------------------------------------------
cross_fitting <- function(
    data,
    K=5,
    proportion = 0.2,
    IV_func = function(data){return(data$Z)},
    print_state = TRUE,
    spline_k = 20,
    constriant = NULL
){
  if(!is.null(constriant)){
    return(cross_fitting_non_deg(data,K,proportion,IV_func,print_state,spline_k))
  }else{
    return(cross_fitting_deg(data,K,proportion,IV_func,print_state,spline_k))
  }
}


cross_fitting_non_deg <- function(
    data,
    K=5,
    proportion = 0.2,
    IV_func = function(data){return(data$Z)},
    print_state = TRUE,
    spline_k = 20
){
  if(is.function(IV_func)){
    # Prespecified an RWF
    delta <- IV_func
    mydat <- data
  }else if(is.numeric(IV_func)){
    # Adaptively train an RWF
    divide = sample(nrow(data), floor(nrow(data)*proportion))
    delta <- delta_adap_func(data[divide,], target_point=IV_func)
    mydat <- data[-divide,]
  }else{
    stop('Error! Please input a valid IV_func!')
  }
  # Define the score waited to be train
  score_AIPW <- rep(0,nrow(mydat))
  score_IPW <- rep(0,nrow(mydat))
  score_OR <- rep(0,nrow(mydat))
  
  score_AIPW_nuc <- rep(0,nrow(mydat))
  score_IPW_nuc <- rep(0,nrow(mydat))
  score_OR_nuc <- rep(0,nrow(mydat))
  
  
  folds <- sample(rep(1:K, length.out = nrow(mydat)))
  for (k in 1:K) {
    train <- which(folds != k)
    divide <- sample(length(train),500)
    train1 <- train[-divide]
    train2 <- train[divide]
    valid <- which(folds == k)
    
    train_set1 <- mydat[train1,]
    train_set2 <- mydat[train2,]
    valid_set <- mydat[valid,]
    
    # Calculate the corresponding delta_ZL
    delta_ZL  <- delta(train_set1)
    
    # Estimate p_{A|L}(a|L)
    L <- mydat$L
    if(length(unique(L))<5){
      pi_model <- npcdensbw(A~factor(L), data = train_set1, bwmethod = "normal-reference")  
    }else{
      pi_model <- npcdensbw(A~L, data = train_set1, bwmethod = "normal-reference")
    }
    pi_m <- npcdens(bws = pi_model)
    
    # Calculate p_{A|L}(a|l) for validation set
    pi_GPS <- predict(pi_m, newdata = valid_set)
    
    # Calculate p_{A}(a) for validation set
    pi_Fix_A_model <- npudensbw(~A, data = train_set1, bwmethod = "normal-reference")
    pi_Fix_A <- fitted(npudens(bws = pi_Fix_A_model, newdata = valid_set))
    
    # Estimate lambda = E[Y|A,L] and rho = E[Z|L]
    if(length(unique(L))<5){
      if(length(unique(valid_set$Y))<=2){
        lambda_model <- gam(Y ~ s(A, by = factor(L), k = spline_k) + factor(L), data=train_set1, 
                            method = "REML", family = binomial())  
      }else{
        lambda_model <- gam(Y ~ s(A, by = factor(L), k = spline_k) + factor(L), data=train_set1, 
                            method = "REML")  
      }  
    }else{
      if (length(unique(valid_set$Y)) <= 2) {
        lambda_model <- gam(Y ~ s(A, L,  k = spline_k), data = train_set1,
                            method = "REML", family = binomial())
      } else {
        lambda_model <- gam(Y ~ s(A, L,  k = spline_k), data = train_set1, 
                            method = "REML")
      }
    }
    lambda <- as.vector(predict(lambda_model, newdata = valid_set, type='response'))
    
    # Estimate E[Z|L]
    if(length(unique(L))<5){
      rho_model <- gam(delta_ZL ~ factor(L), data = train_set1, method = "REML")
    }else{
      rho_model <- gam(delta_ZL ~ s(L, k = spline_k), data = train_set1, method = "REML")
    }
    rho <- as.vector(predict(rho_model, newdata = valid_set, type='response'))
    
    # Estimate E[Z|A,L]
    if(length(unique(L))<5){
      kappa_model <- gam(delta_ZL ~ s(A, by = factor(L), k = spline_k) + factor(L), 
                         data = train_set1, method = "REML")  
    }else{
      kappa_model <- gam(delta_ZL ~ s(A, L,  k = spline_k), 
                         data = train_set1, method = "REML")
    }
    kappa <- as.vector(predict(kappa_model, newdata = valid_set, type='response')) - rho
    kappa <- pmax(pmin(kappa, quantile(kappa,0.95)), quantile(kappa,0.05))
    
    # Estimate mu(A,L) such that E[Y(a)|L]=mu(a,L)
    if(length(unique(L))<5){
      mu_model <- gam(Y * delta_ZL ~ s(A, by = factor(L), k = spline_k) + factor(L), 
                      data = train_set1, method = "REML")  
    }else{
      mu_model <- gam(Y * delta_ZL ~ s(A, L,  k = spline_k), 
                      data = train_set1, method = "REML")
    }
    mu <- as.vector(predict(mu_model, newdata = valid_set, type='response'))
    mu <- (mu - lambda * rho)/kappa
    
    
    
    # Calculate \int mu(A,l) d P_L(l)  for validation set
    grid_data <- do.call(rbind, lapply(valid_set$A, function(a) {
      cbind(A = a, train_set2)
    }))
    
    lambda_vals <- predict(lambda_model, newdata = grid_data, type='response')
    lambda_Fix_A <- tapply(lambda_vals, grid_data$A, mean)
    lambda_Fix_A <- as.numeric(lambda_Fix_A[as.character(valid_set$A)])
    
    
    mu_vals <- predict(mu_model, newdata = grid_data, type = "response") 
    rho_vals <- predict(rho_model, newdata = grid_data, type='response')
    kappa_vals <- predict(kappa_model, newdata = grid_data, type='response') - rho_vals
    kappa_vals <- pmax(pmin(kappa_vals, quantile(kappa_vals,0.95)), quantile(kappa_vals,0.05))
    
    mu_vals <- (mu_vals - lambda_vals * rho_vals)/kappa_vals
    mu_Fix_A <- tapply(mu_vals, grid_data$A, mean)
    mu_Fix_A <- as.numeric(mu_Fix_A[as.character(valid_set$A)])
    
    adj_vals <- (delta(grid_data) - rho_vals) * (lambda_vals-mu_vals)/kappa_vals
    adj_Fix_A <- tapply(adj_vals, grid_data$A, mean)
    adj_Fix_A <- as.numeric(adj_Fix_A[as.character(valid_set$A)])
    
    
    # Calculate the OR, IPW, AIPW score
    delta_ZL_valid <- delta(valid_set)
    score_OR[valid] <- mu_Fix_A
    score_IPW[valid]<- (pi_Fix_A/pi_GPS) * (delta_ZL_valid - rho)* valid_set$Y/kappa
    score_AIPW[valid]<-(pi_Fix_A/pi_GPS) * (delta_ZL_valid - rho)*(valid_set$Y - mu)/kappa + mu_Fix_A - adj_Fix_A
    
    
    score_OR_nuc[valid] <- lambda_Fix_A
    score_IPW_nuc[valid] <- (pi_Fix_A/pi_GPS) * valid_set$Y
    score_AIPW_nuc[valid] <- (pi_Fix_A/pi_GPS) * (valid_set$Y - lambda) + lambda_Fix_A
    
    
    if(print_state == TRUE){
      cat(k,' ')
      if(k==K){cat('\n')}
    }
    
  }
  return(list(
    score = data.frame(
      treatment = mydat$A,
      score_AIPW = score_AIPW,
      score_IPW = score_IPW,
      score_OR = score_OR,
      score_AIPW_nuc = score_AIPW_nuc,
      score_IPW_nuc = score_IPW_nuc,
      score_OR_nuc = score_OR_nuc
    ),
    K=K, spline_k = spline_k, IV_func = IV_func
  ))
}


cross_fitting_deg <- function(
    data,
    K=5, 
    proportion = 0.2,
    IV_func = function(data){return(data$Z)},
    print_state = TRUE,
    spline_k = 20
){
  if(is.function(IV_func)){
    # Prespecified an RWF
    delta <- IV_func
    mydat <- data
  }else if(is.numeric(IV_func)){
    # Adaptively train an RWF
    divide = sample(nrow(data), floor(nrow(data)*proportion))
    delta <- delta_adap_func(data[divide,], target_point=IV_func)
    mydat <- data[-divide,]
  }else{
    stop('Error! Please input a valid IV_func!')
  }
  # Define the score waited to be train
  score_AIPW <- rep(0,nrow(mydat))
  score_IPW <- rep(0,nrow(mydat))
  score_OR <- rep(0,nrow(mydat))
  
  score_AIPW_nuc <- rep(0,nrow(mydat))
  score_IPW_nuc <- rep(0,nrow(mydat))
  score_OR_nuc <- rep(0,nrow(mydat))
  
  
  
  
  folds <- sample(rep(1:K, length.out = nrow(mydat)))
  for (k in 1:K) {
    # k=1
    train <- which(folds != k)
    # divide <- sample(length(train),floor(sqrt(length(train))))
    divide <- sample(length(train),500)
    train1 <- train[-divide]
    train2 <- train[divide]
    valid <- which(folds == k)
    # divide <- sample(length(valid),floor(length(valid)*0.5))
    # valid_1 <- valid[divide]
    # valid_2 <- valid[-divide]
    
    train_set1 <- mydat[train1,]
    train_set2 <- mydat[train2,]
    valid_set <- mydat[valid,]
    
    # Calculate the corresponding delta_ZL
    delta_ZL  <- delta(train_set1)
    
    # Estimate E[Y|A,L]
    Y_type = length(unique(valid_set$Y))
    if(Y_type<=2){
      lambda_model <- gam(Y ~ s(A, k = spline_k), data=train_set1, method = "REML",
                          family = binomial())  
    }else{
      lambda_model <- gam(Y ~ s(A, k = spline_k), data=train_set1, method = "REML")  
    }
    lambda <- as.vector(predict(lambda_model, newdata = valid_set, type='response'))
    
    # Estimate E[Z|L]
    rho_model <- gam(delta_ZL ~ 1, data=train_set1, method = "REML")
    rho <- as.vector(predict(rho_model, newdata = valid_set, type='response'))
    
    # Estimate E[Z|A,L]
    kappa_model <- gam(delta_ZL ~ s(A, k = spline_k), data=train_set1, method = "REML")
    kappa <- as.vector(predict(kappa_model, newdata = valid_set, type='response')) - rho
    kappa <- pmax(pmin(kappa, quantile(kappa,0.95)), quantile(kappa,0.05))
    
    # Estimate mu(A,L) such that E[Y(a)|L]=mu(a,L)
    mu_model <- gam(Y * delta_ZL ~ s(A, k = spline_k), data=train_set1, method = "REML")
    mu <- as.vector(predict(mu_model, newdata = valid_set, type='response'))
    mu <- (mu - lambda * rho)/kappa
    
    
    
    # Calculate \int mu(A,l) d P_L(l)  for validation set
    grid_data <- do.call(rbind, lapply(valid_set$A, function(a) {
      cbind(A = a, train_set2)
    }))
    
    lambda_vals <- predict(lambda_model, newdata = grid_data, type='response')
    lambda_Fix_A <- tapply(lambda_vals, grid_data$A, mean)
    lambda_Fix_A <- as.numeric(lambda_Fix_A[as.character(valid_set$A)])
    
    
    mu_vals <- predict(mu_model, newdata = grid_data, type = "response") 
    rho_vals <- predict(rho_model, newdata = grid_data, type='response')
    kappa_vals <- predict(kappa_model, newdata = grid_data, type='response') - rho_vals
    kappa_vals <- pmax(pmin(kappa_vals, quantile(kappa_vals,0.95)), quantile(kappa_vals,0.05))
    
    mu_vals <- (mu_vals - lambda_vals * rho_vals)/kappa_vals
    mu_Fix_A <- tapply(mu_vals, grid_data$A, mean)
    mu_Fix_A <- as.numeric(mu_Fix_A[as.character(valid_set$A)])
    
    adj_vals <- (delta(grid_data) - rho_vals) * (lambda_vals-mu_vals)/kappa_vals
    adj_Fix_A <- tapply(adj_vals, grid_data$A, mean)
    adj_Fix_A <- as.numeric(adj_Fix_A[as.character(valid_set$A)])
    
    
    # Calculate the OR, IPW, AIPW score
    delta_ZL_valid <- delta(valid_set)
    score_OR[valid] <- mu_Fix_A
    score_IPW[valid]<-  (delta_ZL_valid - rho)* valid_set$Y/kappa
    score_AIPW[valid]<- (delta_ZL_valid - rho)*(valid_set$Y - mu)/kappa + mu_Fix_A - adj_Fix_A
    
    
    score_OR_nuc[valid] <- lambda_Fix_A
    score_IPW_nuc[valid] <- valid_set$Y
    score_AIPW_nuc[valid] <- (valid_set$Y - lambda) + lambda_Fix_A
    
    if(print_state == TRUE){
      cat(k,' ')
      if(k==K){cat('\n')}
    }
    
  }
  
  return(list(
    score = data.frame(
      treatment = mydat$A,
      score_AIPW = score_AIPW,
      score_IPW = score_IPW,
      score_OR = score_OR,
      score_AIPW_nuc = score_AIPW_nuc,
      score_IPW_nuc = score_IPW_nuc,
      score_OR_nuc = score_OR_nuc
    ),
    K=K, spline_k = spline_k, IV_func = IV_func
  ))
}





#' -------------------------------------------------------------------------
#' Function: clean_extreme_scores
#' Purpose: Remove extreme or invalid AIPW scores and optionally restrict to a range of treatment values
#'' Input:
#'   - fit_CF: output of cross_fitting()
#'   - alpha: quantile truncation for score_AIPW (default 0 = no truncation)
#'   - interval: optional range of treatment values c(min, max)
#' Output: list containing cleaned scores and treatment values
#' -------------------------------------------------------------------------
clean_extreme_scores <- function(fit_CF, test, alpha = 0, interval = NULL) {
  # Fixed_A <- fit_CF$IV_test$weak_IV_A
  # Fixed_ind <- fit_CF$IV_test$score <= 1e-3
  Fixed_A <- test$region
  Fixed_ind <- test$pvals <= 0.05
  
  # Determine indices to keep based on weak IV
  ind <- which(pmax(
    approx(Fixed_A, Fixed_ind, xout = fit_CF$score$treatment, method = "constant", rule = 2, f = 0)$y,
    approx(Fixed_A, Fixed_ind, xout = fit_CF$score$treatment, method = "constant", rule = 2, f = 1)$y
  ) == 1)
  
  # Subset scores and treatment
  treatment <- fit_CF$score$treatment[ind]
  score_AIPW <- fit_CF$score$score_AIPW[ind]
  score_IPW <- fit_CF$score$score_IPW[ind]
  score_OR <- fit_CF$score$score_OR[ind]
  score_AIPW_nuc <- fit_CF$score$score_AIPW_nuc[ind]
  score_IPW_nuc <- fit_CF$score$score_IPW_nuc[ind]
  score_OR_nuc <- fit_CF$score$score_OR_nuc[ind]
  
  # Remove NA/Inf
  valid <- is.finite(score_AIPW)
  
  # Quantile truncation
  q_low <- quantile(score_AIPW[valid], alpha, na.rm = TRUE)
  q_high <- quantile(score_AIPW[valid], 1 - alpha, na.rm = TRUE)
  
  if (!is.null(interval)) {
    valid <- valid & score_AIPW >= q_low & score_AIPW <= q_high
    valid <- valid & treatment >= interval[1] & treatment <= interval[2]
  } else {
    valid <- valid & score_AIPW >= q_low & score_AIPW <= q_high
  }
  valid <- valid & treatment >= quantile(Fixed_A,0.05) & treatment <= quantile(Fixed_A,0.95)
  
  return(list(
    valid = valid,
    treatment = treatment[valid],
    score_AIPW = score_AIPW[valid],
    score_IPW = score_IPW[valid],
    score_OR = score_OR[valid],
    score_AIPW_nuc = score_AIPW_nuc[valid],
    score_IPW_nuc = score_IPW_nuc[valid],
    score_OR_nuc = score_OR_nuc[valid]
  ))
}






LLKR <- function(cl){
  #' -------------------------------------------------------------------------
  #' Function: LLKR_func
  #' Purpose: Fit local linear kernel regression using np package
  #' Input:
  #'   - treatment: vector of treatment values
  #'   - score: vector of scores to regress
  #' Output: npreg object representing fitted LLKR
  #' -------------------------------------------------------------------------
  LLKR_func <- function(treatment, score) {
    bw_sub <- npregbw(score ~ treatment, regtype = "ll", ckertype = "epa", bwmethod = "cv.ls", bwscaling = 1)
    LLKR <- npreg(bws = bw_sub, regtype = "ll", ckertype = "epa")
    return(LLKR)
  }
  
  treatment <- cl$treatment
  LLKR_AIPW <- LLKR_func(treatment,cl$score_AIPW)
  LLKR_OR <- LLKR_func(treatment,cl$score_OR)
  LLKR_IPW <- LLKR_func(treatment,cl$score_IPW)
  LLKR_AIPW_nuc <- LLKR_func(treatment,cl$score_AIPW_nuc)
  LLKR_OR_nuc <- LLKR_func(treatment,cl$score_OR_nuc)
  LLKR_IPW_nuc <- LLKR_func(treatment,cl$score_IPW_nuc)
  return(list(
    treatment = treatment,
    LLKR_AIPW = LLKR_AIPW,
    LLKR_IPW = LLKR_IPW,
    LLKR_OR = LLKR_OR,
    LLKR_AIPW_nuc = LLKR_AIPW_nuc,
    LLKR_IPW_nuc = LLKR_IPW_nuc,
    LLKR_OR_nuc = LLKR_OR_nuc
  ))
}

URWF_test <- function(data, func, region = NULL, h = NULL, constriant = NULL, divide = 1){
  if(is.null(constriant)){divide =1}
  if(divide ==1){
    n <- nrow(data); A <- data$A; Z <- func(data)
    if(is.null(h)){h <- sd(A) * n^(-1/4)}
    if(is.null(region)){
      region <- seq(quantile(A,0.05), quantile(A,0.95), length.out = 20)
    }
    if(!is.null(constriant)){
      L <- data[[constriant]]
    }
    res <- numeric(length(region))
    # res2 <- numeric(length(region))
    for(i in 1:length(region)){
      u <- (A - region[i])/h
      Y <- 0.75*(1-u^2)*(abs(u)<=1)/h
      
      if(!is.null(constriant) && length(unique(L))>1){
        if(length(unique(L))<=10){
          fitY <- gam(Y ~ factor(L)); fitZ <- gam(Z ~ factor(L))  
        }else{
          fitY <- gam(Y ~ s(L)); fitZ <- gam(Z ~ s(L))  
        }
        Y <- residuals(fitY); Z <- residuals(fitZ)
      }else{
        Y <- Y - mean(Y); Z <- Z - mean(Z) 
      } 
      # res[i] <- as.numeric(pvalue(independence_test(Y ~ Z, distribution = "approximate")))
      res[i] <- 2*(1-pnorm(abs(mean(Y*Z)/sd(Y*Z) * sqrt(n))))
      # res[i] <- cor.test(Y, Z, method = "spearman")$p.value
      # res2[i] <- cor(Y,Z)
    }
    return(list(region = region, pvals = res))  
  }else{
    n <- nrow(data); A <- data$A; Z <- func(data)
    if(is.null(h)){h <- sd(A) * n^(-1/4)}
    if(is.null(region)){
      region <- seq(quantile(A,0.05), quantile(A,0.95), length.out = 20)
    }
    res <- numeric(length(region))
    data <- data %>% mutate(L_group = ntile(data[[constriant]], divide))
    for(d in 1:divide) {
      res <- URWF_test(data = data %>% filter(L_group==d) %>% select(-L_group), func = func, region = region, 
                       h=h,constriant = constriant, divide = 1)
      if(d==1){
        region <- res$region
        pvals <- res$pvals
      }else{
        pvals <- pmax(pvals, res$pvals)
      }
    }
    return(list(region = region, pvals = pvals))
  }
}

