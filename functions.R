##--------------------------------------------------------------------------------------------------------
##functions for data generation

#X : Unif[0,10] distribution
f_X_unif <- function(n)
{
  return(runif(n,0,10))
}

#Y|X : normal distribution
f_Y_norm <- function(n,mu,sigma)
{
  return(rnorm(n,mu,sigma))
}

#propensity score - setting 1
f_prop_1 <- function(X)
{
  return(1-0.1-0.02*X)
}

#propensity score - setting 2
f_prop_2 <- function(X)
{
  return(1-0.2-0.1*(1+0.1*X)*sin(3*X))
}

##--------------------------------------------------------------------------------------------------------

#function for the computation of weighted quantile
weighted_quantile <- function(probs,values,alpha)
{
  ##inputs
  #-probs : vector of probability masses
  #-values : vector of real numbers
  #-alpha : value in (0,1)
  
  ##output
  #(1-alpha)-th quantile of the discrete distribution
  
  values_sorted <- sort(values)
  q <- min(values_sorted[sapply(values_sorted,function(x) sum(probs[values <= x]))>=1-alpha])
  return(q)
}

##--------------------------------------------------------------------------------------------------------

# Kernel regression function (with Gaussian kernel)
gaussian_kernel <- function(u)
{
  return((1/sqrt(2*pi)) * exp(-(u^2)/2))
}

kernel_regression <- function(X, y, X_pred, bandwidth) {
  n <- nrow(X)
  m <- nrow(X_pred)
  y_pred <- rep(0,m)
  
  for (i in 1:m) 
  {
    weights <- apply(X, 1, function(x_i) gaussian_kernel(sqrt(sum(((x_i - X_pred[i, ])/bandwidth)^2) )))
    y_pred[i] <- sum(weights * y) / sum(weights)
  }
  
  return(y_pred)
}

##--------------------------------------------------------------------------------------------------------
#function for pro-CP
pro_CP <- function(bin,A,score,alpha)
{
  ##inputs
  #-bin : vector of discretized features / bin labels
  #-A : vector of missingness indicators
  #-score : vector of nonconformity scores
  #alpha : target level
  
  ##output
  #(1-alpha)-bound for score, i.e., the (1-alpha)-prediction set C(x) is given by {y : s(x,y) <= output}
  
  bin_unique <- sort(unique(bin))
  n_unique <- length(bin_unique)
  N <- sapply(bin_unique,function(x) sum(bin==x))
  N1 <- sapply(bin_unique,function(x) sum((bin==x)*(A==1)))
  N0 <- sapply(bin_unique,function(x) sum((bin==x)*(A==0)))
  
  bin_observed <- sapply(bin[A==1], function(x) which(bin_unique==x))
  prob_vec <- N0[bin_observed]/(sum(N0)*(N1[bin_observed]+1))
  prob_vec <- c(prob_vec, 1-sum(prob_vec))
  score_vec <- c(score[A==1],1000000)
  return(weighted_quantile(prob_vec,score_vec,alpha))
}

#function for pro-CP2
pro_CP2 <- function(bin,A,score,alpha)
{
  ##inputs
  #-bin : vector of discretized features / bin labels
  #-A : vector of missingness indicators
  #-score : vector of nonconformity scores
  #alpha : target level
  
  ##output
  #(1-alpha)-bound for score, i.e., the (1-alpha)-prediction set C(x) is given by {y : s(x,y) <= output}
  
  bin_unique <- unique(bin)
  M <- length(bin_unique)
  N_vec <- sapply(bin_unique,function(x) sum(bin==x))
  N1_vec <- sapply(bin_unique,function(x) sum((bin==x)*(A==1)))
  N0_vec <- sapply(bin_unique,function(x) sum((bin==x)*(A==0)))
  N0 = sum(N0_vec)
  score_bar <- score
  score_bar[which(A==0)] <- 1000000
  
  probs <- c()
  values <- c()
  for(k in 1:M)
  {
    score_bar_k <- score_bar[which(bin==bin_unique[k])]
    if(N0_vec[k] >= 1)
    {
      probs <- c(probs,rep(N0_vec[k]/N_vec[k],N_vec[k]))
      values <- c(values, score_bar_k)
      
      if(k <= M-1)
      {
        for(j in (k+1):M)
        {
          if(N0_vec[j] >= 1)
          {
            score_bar_j <- score_bar[which(bin==bin_unique[j])]
            probs <- c(probs, rep(2*(N0_vec[k]*N0_vec[j]/(N_vec[k]*N_vec[j])), N_vec[k]*N_vec[j]))
            values <- c(values, c(outer(score_bar_k,score_bar_j,pmin)))
          }
        }
      }
      
    }
    
    if(N0_vec[k] >= 2)
    {
      probs <- c(probs,2*rep(N0_vec[k]*(N0_vec[k]-1)/(N_vec[k]*(N_vec[k]-1)),N_vec[k]*(N_vec[k]-1)/2))
      values <- c(values,apply(combn(1:N_vec[k],2),2,function(v) min(score_bar_k[v])))
    }
    
  }
  
  probs <- probs/(N0^2)
  values_unique  <- sort(unique(score_bar))
  
  probs_simplified <- sapply(values_unique, function(v) sum(probs[values==v]))
  return(weighted_quantile(probs_simplified,values_unique,alpha^2))
  
}

#function for pro-CP with partitioning
pro_CP_U <- function(bin,A,score,alpha,U)
{
  ##inputs
  #-bin : vector of discretized features / bin labels
  #-A : vector of missingness indicators
  #-score : vector of nonconformity scores
  #alpha : target level
  #U : partition
  
  ##output
  #vector of score bounds (length = number of missing labels)
  
  index_observed <- which(A==1)
  index_missing <- which(A==0)
  num_missing <- length(index_missing)
  q_vec <- rep(0,n)
  n_U <- length(U)
  
  for(j in 1:n_U)
  {
    index_missing_j <- intersect(index_missing,U[[j]])
    index_j <- c(index_observed,index_missing_j)
    bin_j <- bin[index_j]
    A_j <- A[index_j]
    score_j <- score[index_j]
    q_j <- pro_CP(bin_j,A_j,score_j,alpha)
    q_vec[index_missing_j] <- rep(q_j,length(index_missing_j))
    
  }
  return(q_vec[index_missing])
}

#function for pro-CP2 with partitioning
pro_CP2_U <- function(bin,A,score,alpha,U)
{
  ##inputs
  #-bin : vector of discretized features / bin labels
  #-A : vector of missingness indicators
  #-score : vector of nonconformity scores
  #alpha : target level
  #U : partition
  
  ##output
  #vector of score bounds (length = number of missing labels)
  
  index_observed <- which(A==1)
  index_missing <- which(A==0)
  num_missing <- length(index_missing)
  q_vec <- c()
  n_U <- length(U)
  
  N0_U <-  sapply(1:n_U, function(j) length(intersect(U[[j]],index_missing)))
  
  for(j in 1:n_U)
  {
    index_missing_j <- intersect(index_missing,U[[j]])
    index_j <- c(index_observed,index_missing_j)
    bin_j <- bin[index_j]
    A_j <- A[index_j]
    score_j <- score[index_j]
    alpha_j <- length(index_missing)*N0_U[j]/(sum(N0_U^2))*alpha
    q_j <- pro_CP2(bin_j,A_j,score_j,alpha)
    q_vec <- c(q_vec,rep(q_j,length(index_missing_j)))
  }
  return(q_vec)
}

##--------------------------------------------------------------------------------------------------------
#function for weighted split conformal prediction
weighted_conformal <- function(A,score,alpha,propensity)
{
  ##inputs
  #-A : vector of missingness indicators
  #-score : vector of nonconformity scores
  #-alpha : target level
  #-propensity : true / estimated propensity score
  
  ##output
  #score bound
  
  m <- sum(1-A)
  w <- (1-propensity)/propensity
  bound <- rep(0,m)
  k <- 1
  for(i in which(A==0))
  {
    indices_i <- c(which(A==1),i)
    weights <- w[indices_i]/sum(w[indices_i])
    values <- c(score[A==1],1000000)
    bound[k] <- weighted_quantile(weights,values,alpha)
    k <- k+1
  }
  return(bound)
}