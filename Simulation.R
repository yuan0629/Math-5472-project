library(mvtnorm)
library(ggplot2)
library(MASS)
library(pracma)
library(patchwork)

# generate data
matrix_generation <- function(row, col, missrate, rank, min, max, sigma) {
  U <- randortho(row, type = "orthonormal")
  V <- randortho(col, type = "orthonormal")
  D <- diag(c(runif(rank, min, max), rep(0, col - rank)), nrow = row, ncol = col)
  X <- U %*% D %*% t(V)
  UV <- X
  SNR <-  sqrt(var(as.vector(X))) / sigma
  E <- matrix(rnorm(col * row, sd = sigma), nrow = row, ncol = col)
  X <- X + E
  miss <- sample(c(1: (col * row)), missrate * col * row, replace = FALSE)
  return(list(X = X, miss = miss, SNR = SNR, UV = UV))
}

# algorithms
soft_impute <- function(X, miss, r, lambda) {
  row <- nrow(X)
  col <- ncol(X)
  Z <- Zold <- Zunobs <- matrix(0, nrow = row, ncol = col)
  RO <- c()
  time <- c()
  
  #t <- 0
  start <- Sys.time()
  while(TRUE){
    #t <- t + 1
    #print(t)
    SVDori <- svd(X+Zunobs)
    U <- SVDori$u
    V <- SVDori$v
    Slambda <- diag(
      sapply(SVDori$d, function(c){max(c-lambda, 0)}),  # hjh: max(vector, 0) not return a vector
      nrow = length(SVDori$d)
    )
    Z <- U %*% Slambda %*% t(V)
    if(sum((Zold - Z)^2) / sum(Z^2) < 1e-6) {
      break
    } else{
      RO <- c(RO, sum((Zold - Z)^2) / sum(Z^2))
      Zold <- Z
      Zunobs <- Z
      Zunobs[-miss] <- 0
      end <- Sys.time()
      time <- c(time, end - start)
    }
  }
  return(list(Z = Z, time = time, RO = RO))
}


SoftImpute_ALS2 <- function(X, miss, r, lambda) {
  row <- nrow(X)
  col <- ncol(X)
  # initialize
  D <- diag(1, nrow = r)
  V <- matrix(0, nrow = col, ncol = r)
  U <- randortho(row, type = "orthonormal")[, 1:r]
  Aold <- A <- U %*% D
  Bold <- B <- V %*% D
  RO <- c()
  time <- c()
  start <- Sys.time()
  
  while(TRUE) {
    
    # update B
    ABunobs <- A %*% t(B)
    ABunobs[-miss] <- 0
    Xstar <- X + ABunobs
    B <- t(Xstar) %*% A %*% solve(t(A) %*% A + lambda * diag(1, nrow = r))
    
    # update A
    ABunobs <- A %*% t(B)
    ABunobs[-miss] <- 0
    Xstar <- X + ABunobs
    A <- Xstar %*% B %*% solve(t(B) %*% B + lambda * diag(1, nrow = r))
    
    
    # convergence
    if(sum((A %*% t(B) - Aold %*% t(Bold))^2) /sum((A %*% t(B))^2) < 1e-6) {
      break
    } else {
      RO <- c(RO, sum((A %*% t(B) - Aold %*% t(Bold))^2) /sum((A %*% t(B))^2))
      Aold <- A
      Bold <- B
      end <- Sys.time()
      time <- c(time, end - start)
    }
  }
  
  return(list(AB = A %*% t(B), RO = RO, time = time))
}

# setting 1
matrix <- matrix_generation(300, 200, 0.7, 50, 1, 500, 5)
X <- matrix$X
miss <- matrix$miss
matrix$SNR
Xobs <- X
Xobs[miss] <- 0

sfals <- SoftImpute_ALS2(Xobs, miss, 25, 120)
timeals <- sfals$time
sum((sfals$AB - X)^2) / sum(X^2)
ROals <- sfals$RO


sf <- soft_impute(Xobs, miss, 25, 120)
sum((sf$Z - X)^2) / sum(X^2)
timesf <- sf$time
ROsf <- sf$RO


p1 <- ggplot() + geom_point(aes(x = timeals, y = ROals, color = "SoftImpute_ALS")) + geom_point(aes(x = timesf, y = ROsf, color = "soft_impute")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) + labs(x = "time", y = "Relative Objective") + 
  theme_minimal() + scale_y_log10() + ggtitle("(300,200) 70%NAs lambda=120 r=15")
p1

# setting 2
matrix <- matrix_generation(800, 600, 0.9, 100, 1, 1000, 5)
X <- matrix$X
miss <- matrix$miss
matrix$SNR
Xobs <- X
Xobs[miss] <- 0

sfals <- SoftImpute_ALS2(Xobs, miss, 50, 140)
timeals <- sfals$time
sum((sfals$AB - X)^2) / sum(X^2)
ROals <- sfals$RO

sf <- soft_impute(Xobs, miss, 50, 140)
sum((sf$Z - X)^2) / sum(X^2)
timesf <- sf$time
ROsf <- sf$RO

p2 <- ggplot() + geom_point(aes(x = timeals, y = ROals, color = "SoftImpute_ALS")) + geom_point(aes(x = timesf, y = ROsf, color = "soft_impute")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) + labs(x = "time", y = "Relative Objective") + 
  theme_minimal() + scale_y_log10() + ggtitle("(800,600) 90%NAs lambda=140 r=50")
p2

# setting 3
matrix <- matrix_generation(1200, 900, 0.8, 100, 1, 1500, 5)
X <- matrix$X
miss <- matrix$miss
matrix$SNR
Xobs <- X
Xobs[miss] <- 0

sfals <- SoftImpute_ALS2(Xobs, miss, 50, 300)
timeals <- sfals$time
sum((sfals$AB - X)^2) / sum(X^2)
ROals <- sfals$RO

sf <- soft_impute(Xobs, miss, 50, 300)
sum((sf$Z - X)^2) / sum(X^2)
timesf <- sf$time
timesf[11:length(timesf)] <- 60 * timesf[11:length(timesf)]
ROsf <- sf$RO

p3 <- ggplot() + geom_point(aes(x = timeals, y = ROals, color = "SoftImpute_ALS")) + geom_point(aes(x = timesf, y = ROsf, color = "soft_impute")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) + labs(x = "time", y = "Relative Objective") + 
  theme_minimal() + scale_y_log10() + ggtitle("(1200,900) 80%NAs lambda=300 r=50")
p3

plot1 <- p1 + p2 + p3 + plot_layout(ncol = 2, byrow = FALSE)
ggsave("D:/personal/homework/Math 5472/project/plot1.png", plot1)

# Simulation2 accuracy

# setting1 SNR=1.7 rank=10
lambda <- seq(from = 1, to = 40, by = 2)
train_errals <- test_errals <- train_errsf <- test_errsf <- rank <- rep(0, length(lambda))
matrix <- matrix_generation(100, 100, 0.5, 10, 1, 100, 1)
X <- matrix$X
miss <- matrix$miss
matrix$SNR
Xobs <- X
Xobs[miss] <- 0
UV <- matrix$UV

for(i in 1:length(lambda)) {
  print(i)
  
  # als
  sfals <- SoftImpute_ALS2(Xobs, miss, 100, lambda[i])
  train_errals[i] <- sum((sfals$AB - X)[-miss]^2) / sum(X[-miss]^2)
  test_errals[i] <- sum((sfals$AB - UV)[miss]^2) / sum(UV[miss]^2)
  
  # sf
  sf <- soft_impute(Xobs, miss, 100, lambda[i])
  train_errsf[i] <- sum((sf$Z - X)[-miss]^2) / sum(X[-miss]^2)
  test_errsf[i] <- sum((sf$Z - UV)[miss]^2) / sum(UV[miss]^2)
  
  rank[i] <- sum(svd(sfals$AB)$d > 1e-2)
}

rank
test_errals
test_errsf

ptest1 <- ggplot() + geom_line(aes(x = rank, y = test_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = test_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "test error") + theme_minimal() + ggtitle("SNR=1.7 rank=10 missing=0.5")

ptrain1 <- ggplot() + geom_line(aes(x = rank, y = train_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = train_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "train error") + theme_minimal()
plot2 <- ptest1 + ptrain1
ggsave("D:/personal/homework/Math 5472/project/plot2.png", plot2, height = 4, width = 9.6)


# setting2 SNR=10 rank=5
lambda <- c(seq(from = 0.1, to = 2, by = 0.2), c(3:12))
train_errals <- test_errals <- train_errsf <- test_errsf <- rank <- rep(0, length(lambda))
matrix <- matrix_generation(100, 100, 0.8, 5, 1, 100, 0.1)
X <- matrix$X
miss <- matrix$miss
matrix$SNR
Xobs <- X
Xobs[miss] <- 0
UV <- matrix$UV
  
for(i in 1:length(lambda)) {
  print(i)
    
  # als
  sfals <- SoftImpute_ALS2(Xobs, miss, 100, lambda[i])
  train_errals[i] <- sum((sfals$AB - X)[-miss]^2) / sum(X[-miss]^2)
  test_errals[i] <- sum((sfals$AB - UV)[miss]^2) / sum(UV[miss]^2)
    
  # sf
  sf <- soft_impute(Xobs, miss, 100, lambda[i])
  train_errsf[i] <- sum((sf$Z - X)[-miss]^2) / sum(X[-miss]^2)
  test_errsf[i] <- sum((sf$Z - UV)[miss]^2) / sum(UV[miss]^2)
    
  rank[i] <- sum(svd(sfals$AB)$d > 1e-2)
}

ptest2 <- ggplot() + geom_line(aes(x = rank, y = test_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = test_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "test error") + theme_minimal() + ggtitle("SNR=17 rank=5 missing=0.8")

ptrain2 <- ggplot() + geom_line(aes(x = rank, y = train_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = train_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "train error") + theme_minimal()
plot3 <- ptest2 + ptrain2
ggsave("D:/personal/homework/Math 5472/project/plot3.png", plot3, height = 4, width = 9.6)


# setting3 SNR=1.7 rank=5
lambda <- seq(from = 1, to = 35, by = 2)
train_errals <- test_errals <- train_errsf <- test_errsf <- rank <- rep(0, length(lambda))
matrix <- matrix_generation(100, 100, 0.5, 5, 1, 100, 1)
X <- matrix$X
miss <- matrix$miss
matrix$SNR
Xobs <- X
Xobs[miss] <- 0
UV <- matrix$UV

for(i in 1:length(lambda)) {
  print(i)
  
  # als
  sfals <- SoftImpute_ALS2(Xobs, miss, 100, lambda[i])
  train_errals[i] <- sum((sfals$AB - X)[-miss]^2) / sum(X[-miss]^2)
  test_errals[i] <- sum((sfals$AB - UV)[miss]^2) / sum(UV[miss]^2)
  
  # sf
  sf <- soft_impute(Xobs, miss, 100, lambda[i])
  train_errsf[i] <- sum((sf$Z - X)[-miss]^2) / sum(X[-miss]^2)
  test_errsf[i] <- sum((sf$Z - UV)[miss]^2) / sum(UV[miss]^2)
  
  rank[i] <- sum(svd(sfals$AB)$d > 1e-2)
}

rank
test_errals
test_errsf

ptest3 <- ggplot() + geom_line(aes(x = rank, y = test_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = test_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "test error") + theme_minimal() + ggtitle("SNR=1.7 rank=5 missing=0.5")

ptrain3 <- ggplot() + geom_line(aes(x = rank, y = train_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = train_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "train error") + theme_minimal()
plot4 <- ptest3 + ptrain3
ggsave("D:/personal/homework/Math 5472/project/plot4.png", plot4, height = 4, width = 9.6)

# setting4 SNR=10 rank=10
lambda <- c(seq(from = 0.1, to = 2, by = 0.2), c(3:12))
train_errals <- test_errals <- train_errsf <- test_errsf <- rank <- rep(0, length(lambda))
matrix <- matrix_generation(100, 100, 0.8, 10, 1, 100, 0.1)
X <- matrix$X
miss <- matrix$miss
matrix$SNR
Xobs <- X
Xobs[miss] <- 0
UV <- matrix$UV

for(i in 1:length(lambda)) {
  print(i)
  
  # als
  sfals <- SoftImpute_ALS2(Xobs, miss, 100, lambda[i])
  train_errals[i] <- sum((sfals$AB - X)[-miss]^2) / sum(X[-miss]^2)
  test_errals[i] <- sum((sfals$AB - UV)[miss]^2) / sum(UV[miss]^2)
  
  # sf
  sf <- soft_impute(Xobs, miss, 100, lambda[i])
  train_errsf[i] <- sum((sf$Z - X)[-miss]^2) / sum(X[-miss]^2)
  test_errsf[i] <- sum((sf$Z - UV)[miss]^2) / sum(UV[miss]^2)
  
  rank[i] <- sum(svd(sfals$AB)$d > 1e-2)
}

ptest4 <- ggplot() + geom_line(aes(x = rank, y = test_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = test_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "test error") + theme_minimal() + ggtitle("SNR=19 rank=10 missing=0.8")

ptrain4 <- ggplot() + geom_line(aes(x = rank, y = train_errals, color = "SoftImpute_ALS"), size = 1) + 
  geom_line(aes(x = rank, y = train_errsf, color = "SoftImpute"), size = 1) + scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  labs(x = "rank", y = "train error") + theme_minimal()
plot5 <- ptest4 + ptrain4
ggsave("D:/personal/homework/Math 5472/project/plot5.png", plot5, height = 4, width = 9.6)



