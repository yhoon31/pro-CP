source("functions.R")

#generate training data
beta0 <- 0
beta1 <- 1
f_X <- f_X_unif
n_train <- 500
X_train <- f_X(n_train)
p_X_train <- 1-0.1-0.02*X_train
#p_X_train <- 1-0.2-0.1*(1+0.1*X_train)*sin(3*X_train)
mu_train = beta0+beta1*X_train
sd_train <- 3+X_train
Y_train <- f_Y_norm(n_train,mu_train,sd_train)
A_train <- rbinom(n_train,1,p_X_train)
betahat <- lm(Y_train[A_train==1]~X_train[A_train==1])$coefficients
#write.table(cbind(X_train,Y_train),"training_data")

##conditional coverage rates of pro-CP and weighted conformal
n <- 500
n_trial <- 500
n_trial_cond <- 100
alpha <- 0.2
eps <- 0.1
K <- 50

#generate partition
U <- list()
for(j in 1:ceiling(n/K))
{
  U[[j]] <- ((j-1)*K+1):min(j*K,n)
}

#run procedures
pro_CP_bin <- matrix(0,nrow=n_trial,ncol=2)
weighted_conformal_bin <- matrix(0,nrow=n_trial,ncol=2)
pro_CP_X <- matrix(0,nrow=n_trial,ncol=2)
weighted_conformal_X <- matrix(0,nrow=n_trial,ncol=2)

pb <- txtProgressBar(min = 0, max = n_trial, style = 3, width = 50, char = "=")

for(t in 1:n_trial)
{
  X <- f_X(n)
  p_X <- 1-0.1-0.02*X
  A <- rbinom(n,1,p_X)
  mu_X = beta0+beta1*X
  sd_X <- 3+X
  muhat <- betahat[1]+betahat[2]*X
  X_discretized_true <- floor(log(p_X/(1-p_X))/log(1+eps))

  pro_CP_bin_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  weighted_conformal_bin_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  pro_CP_X_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  weighted_conformal_X_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  
  for(j in 1:n_trial_cond)
  {
    lower <- (1/0.02)*(0.9-((1+eps)^(X_discretized_true+1))/(1+(1+eps)^(X_discretized_true+1)))
    upper <- (1/0.02)*(0.9-((1+eps)^(X_discretized_true))/(1+(1+eps)^(X_discretized_true)))
    X_trial <- runif(n)
    X_trial <- (upper-lower)*X_trial + lower
    p_X_trial <- 1-0.1-0.02*X_trial
    mu_trial <- beta0+beta1*X_trial
    sd_trial <- 3+X_trial
    Y_bin <- f_Y_norm(n,mu_trial,sd_trial)
    muhat_trial <- betahat[1]+betahat[2]*X_trial
    score_trial <- abs(Y_bin-muhat_trial)
    
    q_1_true <- pro_CP_U(X_discretized_true,A,score_trial,alpha,U)
    q_weighted_conformal_true <- weighted_conformal(A,score_trial,alpha,p_X_trial)
    pro_CP_bin_trial[j,] <- c(mean(score_trial[A==0] <= q_1_true), 2*median(q_1_true))
    weighted_conformal_bin_trial[j,] <- c(mean(score_trial[A==0] <= q_weighted_conformal_true), 2*median(q_weighted_conformal_true))
    
    Y <- f_Y_norm(n,mu_X,sd_X)
    score <- abs(Y-muhat)
    q_1_X <- pro_CP_U(X_discretized_true,A,score,alpha,U)
    q_weighted_conformal_X <- weighted_conformal(A,score,alpha,p_X)
    pro_CP_X_trial[j,] <- c(mean(score[A==0] <= q_1_X), 2*median(q_1_X))
    weighted_conformal_X_trial[j,] <- c(mean(score[A==0] <= q_weighted_conformal_X), 2*median(q_weighted_conformal_X))
    
    
  }
  
  pro_CP_bin[t,] <- c(mean(pro_CP_bin_trial[,1]),mean(pro_CP_bin_trial[,2]))
  weighted_conformal_bin[t,] <- c(mean(weighted_conformal_bin_trial[,1]),mean(weighted_conformal_bin_trial[,2]))
  pro_CP_X[t,] <- c(mean(pro_CP_X_trial[,1]),mean(pro_CP_X_trial[,2]))
  weighted_conformal_X[t,] <- c(mean(weighted_conformal_X_trial[,1]),mean(weighted_conformal_X_trial[,2]))
  
  setTxtProgressBar(pb, t)
}
 

#histograms
cov_breaks <- seq(0.75,0.88,0.003)
len_breaks <- seq(20,28,0.2)

col1 <- rgb(135/255, 206/255, 250/255,alpha=0.5)
col_wc <- rgb(1, 224/255, 197/255,alpha=0.5)

legend1 <- "pro-CP"
legend2 <- "Weighted conformal"

par(mfrow=c(1,2), cex.lab = 1.5, cex.axis=1.5,mar=c(5,5,2,0),oma = c(0, 0, 0, 15),xpd=NA)
hist(pro_CP_bin[,1], breaks = cov_breaks, col = col1, main = "", xlab = "Bin-conditional coverage rate", ylab = "", xlim=c(0.73,0.9), freq=F, ylim=c(0,80))
hist(weighted_conformal_bin[,1], breaks = cov_breaks, col = col_wc, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 80), col = "black", lty = 2, lwd = 2) 

hist(pro_CP_X[,1], breaks = cov_breaks, col = col1, main = "", xlab = "Feature-conditional coverage rate", ylab = "", xlim=c(0.73,0.9), freq=F, ylim=c(0,80))
hist(weighted_conformal_X[,1], breaks = cov_breaks, col = col_wc, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 80), col = "black", lty = 2, lwd = 2)
legend(0.905,45, legend = c(legend1, legend2), fill = c(col1,col_wc),box.lty=0,cex=1.5,bg="transparent", ncol = 1)

##--------------------------------------------------------------------------------------------------------


##feature-conditional coverage and width of pro-CP and weighted conformal

#generate training data
set.seed(1234)
beta0 <- 0
beta1 <- 1
f_X <- f_X_unif
f_prop <- f_prop_1 #setting 1
n_train <- 500
X_train <- f_X(n_train)
p_X_train <- f_prop(X_train)
mu_train = beta0+beta1*X_train
sd_train <- 3+X_train
Y_train <- f_Y_norm(n_train,mu_train,sd_train)
A_train <- rbinom(n_train,1,p_X_train)

#fit regression to construct nonconformity score
betahat <- lm(Y_train[A_train==1]~X_train[A_train==1])$coefficients

#run procedures
n <- 500
n_trial <- 500
n_trial_cond <- 100
alpha <- 0.2
eps <- 0.1

pro_CP_conditional <- matrix(0,nrow=n_trial,ncol=2)
weighted_conformal_conditional <- matrix(0,nrow=n_trial,ncol=2)

pb <- txtProgressBar(min = 0, max = n_trial, style = 3, width = 50, char = "=")

for(t in 1:n_trial)
{
  X <- f_X(n)
  p_X <- f_prop(X)
  A <- rbinom(n,1,p_X)
  mu_t = beta0+beta1*X
  sd_t <- 3+X
  muhat <- betahat[1]+betahat[2]*X
  X_discretized_true <- floor(log(p_X/(1-p_X))/log(1+eps))
  
  pro_CP_conditional_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  weighted_conformal_conditional_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  for(j in 1:n_trial_cond)
  {
    Y <- f_Y_norm(n,mu_t,sd_t)
    score <- abs(Y-muhat)
    q_1_true <- pro_CP_U(X_discretized_true,A,score,alpha,U)
    q_weighted_conformal_true <- weighted_conformal(A,score,alpha,p_X)
    pro_CP_conditional_trial[j,] <- c(mean(score[A==0] <= q_1_true), 2*median(q_1_true))
    weighted_conformal_conditional_trial[j,] <- c(mean(score[A==0] <= q_weighted_conformal_true), 2*median(q_weighted_conformal_true))
    
  }
  
  pro_CP_conditional[t,] <- colMeans(pro_CP_conditional_trial)
  weighted_conformal_conditional[t,] <- colMeans(weighted_conformal_conditional_trial)
  
  setTxtProgressBar(pb, t)
}


#histograms
col1 <- rgb(135/255, 206/255, 250/255,alpha=0.5)
col_wc <- rgb(1, 224/255, 197/255,alpha=0.5)

legend1 <- "pro-CP"
legend2 <- "Weighted conformal"

par(mfrow=c(1,2), cex.lab = 1.5, cex.axis=1.5,mar=c(5,5,2,0),oma = c(0, 0, 0, 15),xpd=NA)
hist(pro_CP_conditional[,1], breaks = seq(0.75,0.9,0.005), col = col1, main = "", xlab = "Bin-conditional coverage rate", ylab = "", xlim=c(0.75,0.9),freq=F, ylim=c(0,50))
hist(weighted_conformal_conditional[,1], breaks = seq(0.75,0.9,0.005), col = col_wc, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 50), col = "black", lty = 2, lwd = 2) 

hist(pro_CP_conditional[,2], breaks = seq(21,27,0.2), col = col1, main = "", xlab = "Feature-conditional coverage rate", ylab = "", freq=F, xlim=c(21,27),ylim=c(0,1.3))
hist(weighted_conformal_conditional[,2], breaks = seq(21,27,0.2), col = col_wc, freq=F, add = TRUE)
legend(27.2,0.7, legend = c(legend1, legend2), fill = c(col1,col_wc),box.lty=0,cex=1.5,bg="transparent", ncol = 1)


