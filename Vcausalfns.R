
# this function computes the probability of each outcome
# from the CPT terms of a V-structure
compute_ps <- function(pX, pZ, pY) {
  ps <- rep(0, 8)
  ps[1] <- (1-pX)*(1-pZ)*(1-pY[1])
  ps[2] <- (1-pX)*(1-pZ)*(pY[1])
  ps[3] <- (1-pX)*(pZ)*(1-pY[2])
  ps[4] <- (1-pX)*(pZ)*(pY[2])
  ps[5] <- (pX)*(1-pZ)*(1-pY[3])
  ps[6] <- (pX)*(1-pZ)*(pY[3])
  ps[7] <- (pX)*(pZ)*(1-pY[4])
  ps[8] <- (pX)*(pZ)*(pY[4])
  return(ps)
}

# add numbers stored in the log space
log_add <- function(es) {
  sum(exp(es - max(es)))*exp(max(es))
}

# add up the series naively for the hypergeometric function
# approx_sd is how many sds to pick either side of binomial mean
naive_3F2 <- function(pX, N, approx_sd = NULL) {
  if (is.null(approx_sd)) {
    range <- 1:N
  } else {
    mn <- pX*N
    sd <- sqrt(pX*(1-pX)*N)
    range_min <- max(round(mn - approx_sd*sd), 1)
    range_max <- min(round(mn + approx_sd*sd), N)
    range <- range_min:range_max
  }
  terms <- lchoose(N, range) - log(range) + range*log(pX) + (N - range)*log(1-pX)
  log_add(terms)
}

# this function computes the theoretical variance of the R estimator
VarR <- function(pX, pZ, pY, N, approx_sd = NULL) {
  ps <- compute_ps(pX, pZ, pY)
  VR1 <- (ps[6] + ps[8])*(ps[5] + ps[7])/(pX^2)*naive_3F2(pX, N, approx_sd)
  VR2 <- (ps[2] + ps[4])*(ps[1] + ps[3])/((1-pX)^2)*naive_3F2(1-pX, N, approx_sd)
  return(VR1 + VR2)
}

# this function computes the theoretical variance of the M estimator
VarM <- function(pX, pZ, pY, N, approx_sd = NULL) {
  ps <- compute_ps(pX, pZ, pY)
  Cterm <- (pY[4] - pY[3] - pY[2] + pY[1])^2*pZ*(1-pZ)/N
  M11factor <- pY[4]*(1-pY[4])
  M11part1 <- M11factor*(ps[3] + ps[4])*naive_3F2(ps[7] + ps[8], N-1, approx_sd)
  M11part2 <- M11factor*(ps[3] + ps[4])^2*(N-1)*naive_3F2(ps[7] + ps[8], N-2, approx_sd)
  M11part3 <- M11factor*(ps[3] + ps[4] + pZ)
  M11 <- (M11part1 + M11part2 + M11part3)/N
  M10factor <- pY[3]*(1-pY[3])
  M10part1 <- M10factor*(ps[1] + ps[2])*naive_3F2(ps[5] + ps[6], N-1, approx_sd)
  M10part2 <- M10factor*(ps[1] + ps[2])^2*(N-1)*naive_3F2(ps[5] + ps[6], N-2, approx_sd)
  M10part3 <- M10factor*(ps[1] + ps[2] + 1 - pZ)
  M10 <- (M10part1 + M10part2 + M10part3)/N
  M01factor <- pY[2]*(1-pY[2])
  M01part1 <- M01factor*(ps[7] + ps[8])*naive_3F2(ps[3] + ps[4], N-1, approx_sd)
  M01part2 <- M01factor*(ps[7] + ps[8])^2*(N-1)*naive_3F2(ps[3] + ps[4], N-2, approx_sd)
  M01part3 <- M01factor*(ps[7] + ps[8] + pZ)
  M01 <- (M01part1 + M01part2 + M01part3)/N
  M00factor <- pY[1]*(1-pY[1])
  M00part1 <- M00factor*(ps[5] + ps[6])*naive_3F2(ps[1] + ps[2], N-1, approx_sd)
  M00part2 <- M00factor*(ps[5] + ps[6])^2*(N-1)*naive_3F2(ps[1] + ps[2], N-2, approx_sd)
  M00part3 <- M00factor*(ps[5] + ps[6] + 1 - pZ)
  M00 <- (M00part1 + M00part2 + M00part3)/N
  return(Cterm + M11 + M10 + M01 + M00)
} 

# this function computes the theoretical covariance of the M and R estimators
CovMR <- function(pX, pZ, pY, N, approx_sd = NULL) {
  ps <- compute_ps(pX, pZ, pY)
  Cov1 <- (ps[5] + ps[7])*(ps[6] + ps[8])/pX^2 + (ps[1] + ps[3])*(ps[2] + ps[4])/(1-pX)^2;
  Fpart1 <- (pY[4]*(1-pY[4])*pZ + pY[3]*(1-pY[3])*(1-pZ))*(1-pX)*naive_3F2(pX, N-1, approx_sd)
  Fpart2 <- (pY[2]*(1-pY[2])*pZ + pY[1]*(1-pY[1])*(1-pZ))*pX*naive_3F2(1-pX, N-1, approx_sd)
  Cov2 <- 2*(pZ - 1)*pZ*(pY[4]-pY[3])*(pY[2]-pY[1]);
  return(Cov1/N + Fpart1 + Fpart2 + Cov2/N)
} 

# this function finds the best linear combination of the M and R estimators
# based on their covariance
find_alpha <- function(pX, pZ, pY, N, approx_sd = NULL) {
  VR <- VarR(pX, pZ, pY, N, approx_sd)
  VM <- VarM(pX, pZ, pY, N, approx_sd)
  CMR <- CovMR(pX, pZ, pY, N, approx_sd)
  alpha <- (VM - CMR) / (VM - CMR + VR - CMR)
  Vmin <- alpha^2*VR + 2*alpha*(1-alpha)*CMR + (1-alpha)^2*VM
  return(list(alpha = alpha, VR = VR, VM = VM, CMR = CMR, Vmin = Vmin))
}

find_alpha_v <- function(VR, VM, CMR) {
  alpha <- (VM - CMR) / (VM - CMR + VR - CMR)
  Vmin <- alpha^2*VR + 2*alpha*(1-alpha)*CMR + (1-alpha)^2*VM
  return(list(alpha = alpha, VR = VR, VM = VM, CMR = CMR, Vmin = Vmin))
}

