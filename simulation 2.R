library(KernSmooth)
source("functions.R")

##--------------------------------------------------------------------------------------------------------

#generate training data / construct score
set.seed(1234)
beta0 <- 0
beta1 <- 1
f_X <- f_X_unif
f_prop <- f_prop_1 #setting 1
#f_prop <- f_prop_2 #setting 2
n_train <- 500
X_train <- f_X(n_train)
p_X_train <- f_prop(X_train) 
mu_train = beta0+beta1*X_train
sd_train <- 3+X_train
Y_train <- f_Y_norm(n_train,mu_train,sd_train)
A_train <- rbinom(n_train,1,p_X_train)
betahat <- lm(Y_train[A_train==1]~X_train[A_train==1])$coefficients ##score = |y - betahat^T x|

#estimate propensity score with kernel regression
h <- dpill(X_train,A_train) ##select bandwidth
phat <- function(x)
{ 
  K_x <- gaussian_kernel((x-X_train)/h)
  return(sum(K_x*A_train)/sum(K_x))
}
phat <- Vectorize(phat)

##--------------------------------------------------------------------------------------------------------
##main simulation
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

sim_est <- matrix(0,nrow=n_trial,ncol=4)
sim_true <- matrix(0,nrow=n_trial,ncol=4)
pb <- txtProgressBar(min = 0, max = n_trial, style = 3, width = 50, char = "=")

for(t in 1:n_trial)
{
  X <- f_X(n)
  p_X <- f_prop(X)
  p_hat <- phat(X)
  
  A <- rbinom(n,1,p_X)
  mu_t = beta0+beta1*X
  sd_t <- 3+X
  Y <- f_Y_norm(n,mu_t,sd_t)
  muhat <- betahat[1]+betahat[2]*X
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

##--------------------------------------------------------------------------------------------------------
col1 <- rgb(135/255, 206/255, 250/255,alpha=0.7)
col2 <- rgb(255/255, 182/255, 193/255,alpha=0.7)
legend1 <- "pro-CP"
legend2 <- "pro-CP2"

par(mfrow=c(2,2),mgp = c(3.5, 1, 0),cex.lab = 1.7, cex.axis=1.7,mar=c(3,5,2,0),oma = c(3, 1, 0, 12),xpd=NA)

hist(sim_true[,1], breaks = seq(0.65,1,0.0075), col = col1, main = "", xlab = "",ylab = "Setting 1", xlim=c(0.65,1), freq=F, ylim=c(0,18))
hist(sim_true[,2], breaks = seq(0.65,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 16), col = "black", lty = 2, lwd = 2) 
hist(sim_true[,3], breaks = seq(18,40,0.5), col = col1, main = "", xlab = "", ylab = "", xlim=c(18,40), freq=F,ylim=c(0,0.25))
hist(sim_true[,4], breaks = seq(18,40,0.5), col = col2, freq=F, add = TRUE)

hist(sim_est[,1], breaks = seq(0.65,1,0.0075), col = col1, main = "", xlab = "", ylab = "Setting 1", xlim=c(0.65,1), freq=F, ylim=c(0,18))
hist(sim_est[,2], breaks = seq(0.65,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 16), col = "black", lty = 2, lwd = 2) 
hist(sim_est[,3], breaks = seq(18,40,0.5), col = col1, main = "", xlab = "", ylab = "", xlim=c(18,40), freq=F,ylim=c(0,0.25))
hist(sim_est[,4], breaks = seq(18,40,0.5), col = col2, freq=F, add = TRUE)

legend(42,0.33, legend = c(legend1, legend2), fill = c(col1,col2),box.lty=0,cex=1.7)

##--------------------------------------------------------------------------------------------------------
##coverage-width curve
set.seed(1234)
f_X <- f_X_unif
n <- 500
n_trial_cond <- 100
eps <- 0.1
K <- 50

alpha_set <- seq(0.15,0.4,by=0.025)
sim1_alpha <- matrix(0,nrow=length(alpha_set),ncol=2)
sim2_alpha <- matrix(0,nrow=length(alpha_set),ncol=2)
k <- 1
  
for(alpha in alpha_set)
{
  sim_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  sim_2_trial <- matrix(0,nrow=n_trial_cond,ncol=2)
  
  for(j in 1:n_trial_cond)
  {
    X_trial <- f_X(n)
    p_X_trial <- f_prop(X_trial)
    A_trial <- rbinom(n,1,p_X_trial)
    X_discretized_trial <- floor(log(p_X_trial/(1-p_X_trial))/log(1+eps))
    mu_trial <- beta0+beta1*X_trial
    sd_trial <- 3+X_trial
    Y_trial <- f_Y_norm(n,mu_trial,sd_trial)
    muhat_trial <- betahat[1]+betahat[2]*X_trial
    score_trial <- abs(Y_trial-muhat_trial)
    
    q_1_true <- pro_CP_U(X_discretized_trial,A_trial,score_trial,alpha,U)
    q_2_true <- pro_CP2_U(X_discretized_trial,A_trial,score_trial,alpha,U)
    sim_trial[j,] <- c(mean(score_trial[A_trial==0] <= q_1_true), 2*median(q_1_true))
    sim_2_trial[j,] <- c(mean(score_trial[A_trial==0] <= q_2_true), 2*median(q_2_true))
    
  } 
  sim1_alpha[k,] <- c(mean(sim_trial[,1]),mean(sim_trial[,2]))
  sim2_alpha[k,] <- c(mean(sim_2_trial[,1]),mean(sim_2_trial[,2]))
  k <- k+1
  
}

par(mfrow=c(1,1),mgp = c(3.5, 1, 0),cex.lab = 1.7, cex.axis=1.7,mar=c(3,5,2,0),oma = c(3, 1, 0, 3))
plot(sim1_alpha,type="l",lty=1,lwd=2,pch=1,xlab="Coverage rate",xlim=c(0.6,0.95),ylim=c(13,33),ylab="Width")
legend("topleft", legend = c("pro-CP", "pro-CP2"), lty=c(1,2),box.lty=0,bg="transparent")
points(sim2_alpha[-1,],type="l",lty=2,lwd=2)

