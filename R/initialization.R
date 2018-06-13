
new_initialization_of_parameters2 <- function(Y, W, V){
  logY <- log(Y)
  W2 <- W^2
  W3 <- W*W2
  W4 <- W*W3
  V2 <- V^2
  V3 <- V*V2
  V4 <- V*V3


  fitlm <- summary(lm(Y ~ V + V2 + V3 + V4 -1))
  Bs <- fitlm$coefficients[, 1]
  B1 <- Bs[1]
  B2 <- Bs[2]
  B3 <- Bs[3]
  B4 <- Bs[4]
  C1 <- B2/B1
  C2 <- B3/B1
  C3 <- B4/B1

  newV <- log(V + C1*V2 + C2*V3 + C3*V4)
  newY <- logY - newV
  as <- summary(lm(newY ~ W + W2 + W3 + W4))$coefficients[, 1]
  a0 <- -0.01
  b1 <- exp(as[1] - a0)
  b2 <- b1*C1
  b3 <- b1*C2
  b4 <- b1*C3
  as[1] <- a0
  bs <- c(b1, b2, b3, b4)
  est <- list(as = as, bs = bs)
  return(est)
}




new_initialization_of_parameters <- function(Y, W, V){
  logY <- log(Y)
  W2 <- W^2
  W3 <- W*W2
  W4 <- W*W3
  V2 <- V^2
  V3 <- V*V2
  V4 <- V*V3

  fitlm <- summary(lm(Y ~ V + V2 + V3 + V4 -1))
  Bs <- fitlm$coefficients[, 1]
  B1 <- Bs[1]
  B2 <- Bs[2]
  B3 <- Bs[3]
  B4 <- Bs[4]
  C1 <- B2/B1
  C2 <- B3/B1
  C3 <- B4/B1
  mV <- mean(V)
  mV2 <- mV*mV
  mV3 <- mV2*mV
  mV4 <- mV3*mV
  F1 <- mV + C1*mV2+ C2*mV3 + C3*mV4

  as <- summary(lm(logY ~ W + W2 + W3 + W4))$coefficients[, 1]
  f1 <- as[1]
  a1 <- as[2]
  a2 <- as[3]
  a3 <- as[4]
  a4 <- as[5]
  mW <- mean(W)
  mW2 <- mW*mW
  mW3 <- mW2*mW
  mW4 <- mW3*mW
  F2 <- a1*mW + a2*mW2 + a3*mW3 + a4*mW4
  a0 <- (log(B1) + f1 - F2 - log(F1))/2
  b1 <- exp(f1 - a0)/F1
  b2 <- b1*C1
  b3 <- b1*C2
  b4 <- b1*C3

  newV <- log(b1*V + b2*V2 + b3*V3 + b4*V4)
  newY <- logY - newV
  as <- summary(lm(newY ~ W + W2 + W3 + W4))$coefficients[, 1]

  bs <- c(b1, b2, b3, b4)
  est <- list(as = as, bs = bs)
  return(est)
}




correct_individual_sample <- function(Y, gc, reads, gamma =  0.8, steps = 50, down = 0.1){

  data <- Y
  flag <- Y!=0
  Y <- Y[flag]
  gc <- gc[flag]
  reads <- reads[flag]
  est <- new_initialization_of_parameters(Y, gc, reads)
  as <- unlist(est$as)
  bs <- unlist(est$bs)
  a0 <- as[1]
  a1 <- as[2]
  a2 <- as[3]
  a3 <- as[4]
  a4 <- as[5]
  b1 <- bs[1]
  b2 <- bs[2]
  b3 <- bs[3]
  b4 <- bs[4]
  res <- gradient_descent_for_individual_sample(Y, gc, reads, a0, a1, a2, a3, a4, b1, b2, b3, b4, gamma, steps, down)

  data[flag] <- res$corrected
  report <- list(parameter = res$parameter, corrected = data)

  return(report)
}
