# This code simulates binary data from a V-structure
# and evaluates different causal effect estimators

pX <- 1/3 # probability of node X
pZ <- 2/3 # probability of node Z
pY <- c(1/6, 1/2, 1/3, 5/6) # probabilities of node Y given X and Z

N <- 100 # sample size
nreps <- 4e7 # number of MC iterations

set.seed(42)

# probabilities of the binary vectors
p0 <- (1-pX)*(1-pZ)*(1-pY[1])
p1 <- (1-pX)*(1-pZ)*(pY[1])
p2 <- (1-pX)*(pZ)*(1-pY[2])
p3 <- (1-pX)*(pZ)*(pY[2])
p4 <- (pX)*(1-pZ)*(1-pY[3])
p5 <- (pX)*(1-pZ)*(pY[3])
p6 <- (pX)*(pZ)*(1-pY[4])
p7 <- (pX)*(pZ)*(pY[4])

R <- rep(NA, nreps)
M <- rep(NA, nreps)

for (ii in 1:nreps) {
  # sample binary data, we shift the labelling by 1
  binary_sample <- sample(1:8, N, replace = TRUE, prob = c(p0, p1, p2, p3, p4, p5, p6, p7))
  # indexing is still shifted by 1
  Ns <- tabulate(binary_sample, 8)
  # compute R and M estimators accounting for index shift
  Rparts <- c((Ns[2]+Ns[4])/sum(Ns[1:4]), (Ns[6]+Ns[8])/sum(Ns[5:8]))
  # set empty count parts to 0
  Rparts[which(is.na(Rparts))] <- 0
  R[ii] <- Rparts[2] - Rparts[1]
  Mparts <- c(Ns[2]/sum(Ns[1:2]), Ns[4]/sum(Ns[3:4]), 
              Ns[6]/sum(Ns[5:6]), Ns[8]/sum(Ns[7:8]))
  # set empty count parts to 0
  Mparts[which(is.na(Mparts))] <- 0
  M[ii] <- (Mparts[4] - Mparts[2])*(sum(Ns[3:4]) + sum(Ns[7:8]))/N +
           (Mparts[3] - Mparts[1])*(sum(Ns[1:2]) + sum(Ns[5:6]))/N
}

# causal effect
print(mean(R))
print(mean(M))

# sd of causal estimators
print(sd(R))
print(sd(M))
print(sqrt(cov(R,M)))

# load R functions to compute theoretical variances
source("Vcausalfns.R")
# theoretical sd of causal estimators 
opt <- find_alpha(pX, pZ, pY, N)
print(sqrt(opt$VR))
print(sqrt(opt$VM))
print(sqrt(opt$CMR))

# optimal combination parameter
print(opt$alpha)
# causal effect from optimal combination
print(mean(opt$alpha*R + (1-opt$alpha)*M))
# sd of combination
print(sd(opt$alpha*R + (1-opt$alpha)*M))
# theoretical optimum
print(sqrt(opt$Vmin))
