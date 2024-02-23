source("functions.R")
library(tidyverse)

##--------------------------------------------------------------------------------------------------------
#load data
data <- get(load("Jobs-NoMiss-cont.RData"))
data_treated <- data[which(data$treat=="exp"),]
dim(data_treated)

X <- data_treated[,-c(1,12,14)]
Y <- data_treated$depress2

#data preprocessing
X = mutate(X, new_col = rep(1,dim(X)[1]))
X <- pivot_wider(X, names_from = c("occp"), values_from="new_col")
X = mutate(X, new_col = rep(1,dim(X)[1]))
X <- pivot_wider(X, names_from = c("marital"), values_from="new_col")
X <- mutate(X, educ = recode(educ, "lt-hs"=0, "highsc"=1, "somcol"=2, "bach" = 3, "gradwk" = 4))
X <- mutate(X, nonwhite = recode(nonwhite, "white0"=0, "non.white1"=1))
X <- mutate(X, income = recode(income, "lt15k"=0, "15t24k"=1, "25t39k"=2, "40t49k" = 3, "50k+" = 4))
X <- mutate(X, job_disc = recode(job_disc, "Low-Med"=0, "Medium"=1, "Med-High"=2, "High" = 3))
X <- X %>% mutate_all(~ifelse(is.na(.), 0, .))
X <- as.matrix(X)

##--------------------------------------------------------------------------------------------------------
#generate train & calibration data with data splitting
set.seed(1234)
n <- dim(X)[1]
n_train <- 379
n_cal <- n-n_train #500
ind_train <- sort(sample(1:n,n_train,replace=F))
ind_cal <- setdiff(1:n,ind_train)
X_train <- X[ind_train,]
Y_train <- Y[ind_train]
X_cal <- X[ind_cal,]
Y_cal <- Y[ind_cal]

#simulate missingness
beta_true <- c(-0.05,0.05,-0.05,0.03,0.03,-0.1,-0.1,-0.1,0.03,0.05,0.03,0.05,-0.1,-0.1,0.05,-0.1,0.1,-0.1,-0.1,-0.1,-0.2,-0.1,-0.2)
beta0_true <- 1
p_A_train <- exp(beta0_true+X_train%*%beta_true)/(1+exp(beta0_true+X_train%*%beta_true))

mean(p_A_train)
max(p_A_train)
A_train <- rbinom(n_train,1,p_A_train)
plot(sort(p_A_train),ylim=c(0.4,1))

p_A <- exp(beta0_true+X_cal%*%beta_true)/(1+exp(beta0_true+X_cal%*%beta_true))

#construct and compute nonconformity score
library(randomForest)
rf <- randomForest(X_train[which(A_train==1),],Y_train[which(A_train==1)])
muhat <- unname(predict(rf,X_cal))
score <- abs(Y_cal - muhat)

#estimate propensity score
rf_propensity <- randomForest(X_train,as.factor(A_train))
p_hat <- predict(rf_propensity,X_cal,type="prob")[,2]

##--------------------------------------------------------------------------------------------------------
#run procedures
alpha <- 0.2
eps <- 0.1
n_trial <- 500

X_discretized_true <- floor(log(p_A/(1-p_A))/log(1+eps))
X_discretized_est <- floor(log(p_hat/(1-p_hat))/log(1+eps))

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
  q_1 <- pro_CP_U(X_discretized_est,A,score_observed,alpha,U)
  q_2 <- pro_CP2_U(X_discretized_est,A,score_observed,alpha,U)
  
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

##--------------------------------------------------------------------------------------------------------
#histograms

col1 <- rgb(135/255, 206/255, 250/255,alpha=0.7)
col2 <- rgb(255/255, 182/255, 193/255,alpha=0.7)
legend1 <- "pro-CP"
legend2 <- "pro-CP2"

par(mfrow=c(2,2),mgp = c(3.5, 1, 0),cex.lab = 1.7, cex.axis=1.7,mar=c(3,5,2,0),oma = c(3, 1, 0, 12),xpd=NA)

hist(coverage_true_1, breaks = seq(0.67,1,0.0075), col = col1, main = "", xlab = "", ylab = "Known propensity score", xlim=c(0.67,1), freq=F, ylim=c(0,15))
hist(coverage_true_2, breaks = seq(0.67,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 15), col = "black", lty = 2, lwd = 2) 
hist(width_true_1, breaks = seq(1.1,2.35,0.02), col = col1, main = "", xlab = "", ylab = "", xlim=c(1.1,2.35), freq=F,ylim=c(0,10))
hist(width_true_2, breaks = seq(1.1,2.35,0.02), col = col2, freq=F, add = TRUE)

hist(coverage_1, breaks = seq(0.67,1,0.0075), col = col1, main = "", xlab = "Coverage proportion", ylab = "Unknown propensity score", xlim=c(0.67,1), freq=F, ylim=c(0,15))
hist(coverage_2, breaks = seq(0.67,1,0.0075), col = col2, freq=F, add = TRUE)
lines(x = c(0.8, 0.8), y = c(0, 15), col = "black", lty = 2, lwd = 2) 
hist(width_1, breaks = seq(1.1,2.7,0.03), col = col1, main = "", xlab = "Median width", ylab = "", xlim=c(1.1,2.7), freq=F,ylim=c(0,8))
hist(width_2, breaks = seq(1.1,2.7,0.03), col = col2, freq=F, add = TRUE)

legend(2.85,11, legend = c(legend1, legend2), fill = c(col1,col2),box.lty=0,cex=1.7)




