# this function computes the probability of each outcome
# from the CPT terms of a V-structure with fixed X
compute_qs <- function(pZ, pY) {
  qs <- rep(0, 8)
  qs[1] <- (1-pZ)*(1-pY[1])
  qs[2] <- (1-pZ)*(pY[1])
  qs[3] <- (pZ)*(1-pY[2])
  qs[4] <- (pZ)*(pY[2])
  qs[5] <- (1-pZ)*(1-pY[3])
  qs[6] <- (1-pZ)*(pY[3])
  qs[7] <- (pZ)*(1-pY[4])
  qs[8] <- (pZ)*(pY[4])
  return(qs)
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

# this function computes the theoretical variance of the R estimator with fixed X
VarR_fixedX <- function(pZ, pY, Ns) {
  qs <- compute_qs(pZ, pY)
  VR1 <- (qs[6] + qs[8])*(qs[5] + qs[7])/Ns[2]
  VR2 <- (qs[2] + qs[4])*(qs[1] + qs[3])/Ns[1]
  return(VR1 + VR2)
}

# this function computes the theoretical variance of the M estimator with fixed X
VarM_fixedX <- function(pZ, pY, Ns, approx_sd = NULL) {
  qs <- compute_qs(pZ, pY)
  N <- sum(Ns)
  Cterm11 <- N*qs[8]*(1-qs[8]) + Ns[1]*qs[7]*qs[8]/pZ
  M11 = Ns[1]*qs[7]*qs[8]*(1 + (Ns[1]-1)*pZ)/(pZ)*naive_3F2(pZ, Ns[2], approx_sd)
  Cterm10 <- N*qs[6]*(1-qs[6]) + Ns[1]*qs[5]*qs[6]/(1-pZ)
  M10 = Ns[1]*qs[5]*qs[6]*(1 + (Ns[1]-1)*(1-pZ))/(1-pZ)*naive_3F2(1-pZ, Ns[2], approx_sd)
  Cterm01 <- N*qs[4]*(1-qs[4]) + Ns[2]*qs[3]*qs[4]/pZ
  M01 = Ns[2]*qs[3]*qs[4]*(1 + (Ns[2]-1)*pZ)/(pZ)*naive_3F2(pZ, Ns[1], approx_sd)
  Cterm00 <- N*qs[2]*(1-qs[2]) + Ns[2]*qs[1]*qs[2]/(1-pZ)
  M00 = Ns[2]*qs[1]*qs[2]*(1 + (Ns[2]-1)*(1-pZ))/(1-pZ)*naive_3F2(1-pZ, Ns[1], approx_sd)
  Covpart1 = (qs[8] - qs[4])*(qs[2] - qs[6])
  Covpart2 = qs[8]*qs[4]*(1-pZ)/pZ + qs[2]*qs[6]*pZ/(1-pZ)
  Varpart = (Cterm11 + M11 + Cterm10 + M10 + Cterm01 + M01 + Cterm00 + M00)/N^2
  Covpart = 2*(Covpart1-Covpart2)/N
  return(Varpart + Covpart)
} 

# this function computes the theoretical covariance of the M and R estimator with fixed X
CovMR_fixedX <- function(pZ, pY, Ns) {
  qs <- compute_qs(pZ, pY)
  N <- sum(Ns)
  part1 <- -(qs[8] + qs[6] - qs[4] - qs[2])^2/N
  part2 <- (qs[8] - qs[4])^2/(pZ*N) + (qs[6] - qs[2])^2/((1-pZ)*N)
  part3 <- qs[7]*qs[8]/(pZ*Ns[2]) + qs[5]*qs[6]/((1-pZ)*Ns[2]) + qs[3]*qs[4]/(pZ*Ns[1]) + qs[1]*qs[2]/((1-pZ)*Ns[1])
  cov_result <- part1 + part2 + part3
  return(cov_result)
}

# this function finds the best linear combination of the M and R estimators
# based on their covariance
find_alpha_fixedX <- function(pZ, pY, Ns, approx_sd = NULL) {
  VR <- VarR_fixedX(pZ, pY, Ns)
  VM <- VarM_fixedX(pZ, pY, Ns, approx_sd)
  CMR <- CovMR_fixedX(pZ, pY, Ns)
  alpha <- (VM - CMR) / (VM - CMR + VR - CMR)
  Vmin <- alpha^2*VR + 2*alpha*(1-alpha)*CMR + (1-alpha)^2*VM
  return(list(alpha = alpha, VR = VR, VM = VM, CMR = CMR, Vmin = Vmin))
}

find_alpha_v <- function(VR, VM, CMR) {
  alpha <- (VM - CMR) / (VM - CMR + VR - CMR)
  Vmin <- alpha^2*VR + 2*alpha*(1-alpha)*CMR + (1-alpha)^2*VM
  return(list(alpha = alpha, VR = VR, VM = VM, CMR = CMR, Vmin = Vmin))
}