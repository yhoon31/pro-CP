library(randomForest)
library(MASS)
source("functions.R")

#generate training data
set.seed(1234)
p <- 30
mu <- rep(1,p)
sigma <- diag(2,p)
beta0 <- 5
beta <- runif(p,-2,2)
n_train <- 500
X_train <- mvrnorm(n_train,mu,sigma)
beta_propensity <- runif(p,-0.1,0.1)
beta0_propensity <- 1.5
p_X_train <- exp(beta0_propensity+X_train%*%beta_propensity)/(1+exp(beta0_propensity+X_train%*%beta_propensity)) #logistic model
mean(p_X_train) #~24% missingness
mu_train = beta0+X_train%*%beta
sd_train <- rowMeans(X_train^2)
Y_train <- f_Y_norm(n_train,mu_train,sd_train)
A_train <- rbinom(n_train,1,p_X_train)
betahat <- lm(Y_train[A_train==1]~X_train[A_train==1,])$coefficients #fit regression to construct score

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

#main simulation
pro_CP_bin <- matrix(0,nrow=n_trial,ncol=2)
weighted_conformal_bin <- matrix(0,nrow=n_trial,ncol=2)
pro_CP_X <- matrix(0,nrow=n_trial,ncol=2)
weighted_conformal_X <- matrix(0,nrow=n_trial,ncol=2)

pb <- txtProgressBar(min = 0, max = n_trial, style = 3, width = 50, char = "=")

for(t in 1:n_trial)
{
  X <- mvrnorm(n,mu,sigma)
  p_X <- exp(beta0_propensity+X%*%beta_propensity)/(1+exp(beta0_propensity+X%*%beta_propensity))
  A <- rbinom(n,1,p_X)
  mu_X = beta0+X%*%beta
  sd_X <- rowMeans(X^2)
  muhat <- betahat[1]+X%*%betahat[-1]
  X_discretized_true <- floor(log(p_X/(1-p_X))/log(1+eps))
  
  pro_CP_X_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  weighted_conformal_X_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  
  for(j in 1:n_trial_cond)
  {
    Y <- f_Y_norm(n,mu_X,sd_X)
    score <- abs(Y-muhat)
    q_1_X <- pro_CP_U(X_discretized_true,A,score,alpha,U)
    q_weighted_conformal_X <- weighted_conformal(A,score,alpha,p_X)
    pro_CP_X_trial[j,] <- c(mean(score[A==0] <= q_1_X), 2*median(q_1_X))
    weighted_conformal_X_trial[j,] <- c(mean(score[A==0] <= q_weighted_conformal_X), 2*median(q_weighted_conformal_X))
  }
  
  pro_CP_X[t,] <- c(mean(pro_CP_X_trial[,1]),mean(pro_CP_X_trial[,2]))
  weighted_conformal_X[t,] <- c(mean(weighted_conformal_X_trial[,1]),mean(weighted_conformal_X_trial[,2]))
  
  setTxtProgressBar(pb, t)
}

#plots
col1 <- rgb(135/255, 206/255, 250/255,alpha=0.7)
col_wc <- rgb(1, 224/255, 197/255,alpha=0.7)
legend1 <- "pro-CP"
legend2 <- "weighted conformal"

par(mfrow=c(1,2),mgp = c(3.5, 1, 0),cex.lab = 1.7, cex.axis=1.7,mar=c(3,5,2,0),oma = c(3, 0, 0, 17),xpd=NA)

hist(pro_CP_X[,1], breaks = seq(0.75,1,0.005), col = col1, main = "", xlab = "Feature-conditional coverage rate", ylab="", xlim=c(0.75,1), freq=F, ylim=c(0,35))
hist(weighted_conformal_X[,1], breaks = seq(0.75,1,0.005), col = col_wc, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 35), col = "black", lty = 2, lwd = 2) 
hist(pro_CP_X[,2], breaks = seq(7.5,12.5,0.1), col = col1, main = "", xlab = "Conditional expectation of median width", ylab = "", xlim=c(7.5,12.5), freq=F,ylim=c(0,4))
hist(weighted_conformal_X[,2], breaks = seq(7.5,12.5,0.1), col = col_wc, freq=F, add = TRUE)

legend(12.6,2, legend = c(legend1, legend2), fill = c(col1,col_wc),box.lty=0,cex=1.7)


##illustrate performance of pro-CP and pro-CP2

#estimate propensity score
rf <- randomForest(X_train,as.factor(A_train))

set.seed(1234)
n <- 500
n_trial <- 500
alpha <- 0.2
eps <- 0.1

#construct partition
K <- 50

U <- list()
for(j in 1:ceiling(n/K))
{
  U[[j]] <- ((j-1)*K+1):min(j*K,n)
}

#main simulation
sim_est <- matrix(0,nrow=n_trial,ncol=4)
sim_true <- matrix(0,nrow=n_trial,ncol=4)
pb <- txtProgressBar(min = 0, max = n_trial, style = 3, width = 50, char = "=")

for(t in 1:n_trial)
{
  X <- mvrnorm(n,mu,sigma)
  p_X <- exp(beta0_propensity+X%*%beta_propensity)/(1+exp(beta0_propensity+X%*%beta_propensity))
  p_hat <- predict(rf,X,type="prob")[,2]
  
  A <- rbinom(n,1,p_X)
  mu_X = beta0+X%*%beta
  sd_X <- rowMeans(X^2)
  muhat <- betahat[1]+X%*%betahat[-1]
  Y <- f_Y_norm(n,mu_X,sd_X)
  score <- abs(Y-muhat)
  
  X_discretized <- floor(log(p_hat/(1-p_hat))/log(1+eps))
  
  q_1 <- pro_CP_U(X_discretized,A,score,alpha,U)
  q_2 <- pro_CP2_U(X_discretized,A,score,alpha,U)
  sim_est[t,] <- c(mean(score[A==0] <= q_1), mean(score[A==0] <= q_2), 2*median(q_1), 2*median(q_2))
  
  X_discretized_true <- floor(log(p_X/(1-p_X))/log(1+eps))
  q_1_true <- pro_CP_U(X_discretized_true,A,score,alpha,U)
  q_2_true <- pro_CP2_U(X_discretized_true,A,score,alpha,U)
  sim_true[t,] <- c(mean(score[A==0] <= q_1_true), mean(score[A==0] <= q_2_true), 2*median(q_1_true), 2*median(q_2_true))
  
  
  setTxtProgressBar(pb, t)
}

#plots
col1 <- rgb(135/255, 206/255, 250/255,alpha=0.7)
col2 <- rgb(255/255, 182/255, 193/255,alpha=0.7)
legend1 <- "pro-CP"
legend2 <- "pro-CP2"

par(mfrow=c(2,2),mgp = c(3.5, 1, 0),cex.lab = 1.7, cex.axis=1.7,mar=c(3,5,2,0),oma = c(3, 1, 0, 12),xpd=NA)

hist(sim_true[,1], breaks = seq(0.7,1,0.0075), col = col1, main = "", xlab = "", ylab = "Known propensity score", xlim=c(0.7,1), freq=F, ylim=c(0,14))
hist(sim_true[,2], breaks = seq(0.7,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 14), col = "black", lty = 2, lwd = 2) 
hist(sim_true[,3], breaks = seq(7.5,16,0.2), col = col1, main = "", xlab = "", ylab = "", xlim=c(7.5,16), freq=F,ylim=c(0,0.8))
hist(sim_true[,4], breaks = seq(7.5,16,0.2), col = col2, freq=F, add = TRUE)
hist(sim_est[,1], breaks = seq(0.7,1,0.0075), col = col1, main = "", xlab = "Coverage proportion", ylab = "Unknown propensity score", xlim=c(0.7,1), freq=F, ylim=c(0,14))
hist(sim_est[,2], breaks = seq(0.7,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 14), col = "black", lty = 2, lwd = 2) 
hist(sim_est[,3], breaks = seq(7.5,16,0.2), col = col1, main = "", xlab = "Median width", ylab = "", xlim=c(7.5,16), freq=F,ylim=c(0,0.7))
hist(sim_est[,4], breaks = seq(7.5,16,0.2), col = col2, freq=F, add = TRUE)

legend(16.75,0.9, legend = c(legend1, legend2), fill = c(col1,col2),box.lty=0,cex=1.7)


