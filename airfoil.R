library(randomForest)
source("functions.R")

# load data
data = read.table("airfoil.txt")
colnames(data) = c("Frequency","Angle","Chord","Velocity","Suction","Sound")

X = as.matrix(data[,1:5])
Y = as.numeric(data[,6])
X[,1] = log(X[,1]) # Log transform
X[,5] = log(X[,5]) # Log transform
n = nrow(X)
p = ncol(X)

#data splitting
set.seed(1234)
n_train <- 803
n_cal <- n-n_train #700
ind_train <- sort(sample(1:n,n_train,replace=F))
ind_cal <- setdiff(1:n,ind_train)
X_train <- X[ind_train,]
Y_train <- Y[ind_train]
X_cal <- X[ind_cal,]
Y_cal <- Y[ind_cal]


##simulate missingness
set.seed(1234)
beta_true <- c(0.1,0.05,1,0.01,0.1)
p_A_train <- exp(X_train%*%beta_true)/(1+exp(X_train%*%beta_true))
A_train <- rbinom(n_train,1,p_A_train)

##construct nonconformity score with random forest regression
rf <- randomForest(X_train[which(A_train==1),],Y_train[which(A_train==1)])
muhat <- unname(predict(rf,X_cal))
score <- abs(Y_cal - muhat)


#true propensity score
p_A <- exp(X_cal%*%beta_true)/(1+exp(X_cal%*%beta_true))

##estimate propensity score
set.seed(1234)
#select bandwidth with single validation
index_1 <- sample(dim(X_train)[1],400,replace=F)
X_train_1 <- X_train[index_1,]
X_train_2 <- X_train[-index_1,]
A_train_1 <- A_train[index_1]
A_train_2 <- A_train[-index_1]
X_sd <- apply(X_train_1,2,sd)
h_set <- seq(1,3,0.25)
err <- rep(0,length(h_set))
l <- 1
for(h in h_set)
{
  p_hat_h <- kernel_regression(X_train_1,A_train_1,X_train_2,X_sd*h)
  err[l] <- sqrt(sum((A_train_2 - p_hat_h)^2))
  l <- l+1
}
bandwidth <- X_sd*h_set[which.min(err)]
p_hat <- kernel_regression(X_train,A_train,X_cal,bandwidth)


##run procedures
alpha <- 0.2
eps <- 0.1
n_trial <- 500
X_discretized_true <- floor(log(p_A/(1-p_A))/log(1+eps))
X_discretized <- floor(log(p_hat/(1-p_hat))/log(1+eps))

#construct partition
K <- 50 #partition size 14
U <- list()
for(j in 1:ceiling(n_cal/K))
{
  U[[j]] <- ((j-1)*K+1):min(j*K,n_cal)
}

coverage_true_1 <- rep(0,n_trial)
coverage_true_2 <- rep(0,n_trial)
width_true_1 <- rep(0,n_trial)
width_true_2 <- rep(0,n_trial)
coverage_1 <- rep(0,n_trial)
coverage_2 <- rep(0,n_trial)
width_1 <- rep(0,n_trial)
width_2 <- rep(0,n_trial)

pb <- txtProgressBar(min = 0, max = n_trial, style = 3,  width = 50, char = "=")
for(i in 1:n_trial)
{
  A <- rbinom(n_cal,1,p_A)
  score_observed <- score*A
  q_1_true <- pro_CP_U(X_discretized_true,A,score_observed,alpha,U)
  q_2_true <- pro_CP2_U(X_discretized_true,A,score_observed,alpha,U)
  q_1 <- pro_CP_U(X_discretized,A,score_observed,alpha,U)
  q_2 <- pro_CP2_U(X_discretized,A,score_observed,alpha,U)
  
  score_missing <- score[which(A==0)]
  
  coverage_true_1[i] <- mean(score_missing <= q_1_true)
  width_true_1[i] <- 2*median(q_1_true)
  coverage_true_2[i] <- mean(score_missing <= q_2_true)
  width_true_2[i] <- 2*median(q_2_true)
  coverage_1[i] <- mean(score_missing <= q_1)
  width_1[i] <- 2*median(q_1)
  coverage_2[i] <- mean(score_missing <= q_2)
  width_2[i] <- 2*median(q_2)
  
  setTxtProgressBar(pb, i)
}

#histograms
col1 <- rgb(135/255, 206/255, 250/255,alpha=0.7)
col2 <- rgb(255/255, 182/255, 193/255,alpha=0.7)
legend1 <- "pro-CP"
legend2 <- "pro-CP2"

par(mfrow=c(2,2),cex.lab = 1.6, cex.axis=1.5,mar=c(5,5,2,2))

hist(coverage_true_1, breaks = seq(0.7,1,0.0075), col = col1, main = "", xlab = "", ylab = "Known propensity score", xlim=c(0.7,1), freq=F, ylim=c(0,20))
hist(coverage_true_2, breaks = seq(0.7,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 20), col = "black", lty = 2, lwd = 2) 
hist(width_true_1, breaks = seq(8,13,0.12), col = col1, main = "", xlab = "", ylab = "", xlim=c(8,13), freq=F,ylim=c(0,3))
hist(width_true_2, breaks = seq(8,13,0.12), col = col2, freq=F, add = TRUE)

hist(coverage_1, breaks = seq(0.7,1,0.0075), col = col1, main = "", xlab = "Coverage proportion", ylab = "Unknown propensity score", xlim=c(0.7,1), freq=F, ylim=c(0,20))
hist(coverage_2, breaks = seq(0.7,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 20), col = "black", lty = 2, lwd = 2) 
hist(width_1, breaks = seq(8,13,0.12), col = col1, main = "", xlab = "Median width", ylab = "", xlim=c(8,13), freq=F,ylim=c(0,4))
hist(width_2, breaks = seq(8,13,0.12), col = col2, freq=F, add = TRUE)

legend(13.5,5.8, legend = c(legend1, legend2), fill = c(col1,col2),box.lty=0,cex=1.7)
