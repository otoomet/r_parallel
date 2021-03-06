
library(pbdMPI, quiet=TRUE)
init()
rank <- comm.rank()
master <- comm.rank() == 0
comm.cat("Job start", as.character(date()), "\n", quiet=TRUE)
##


mu <- 1
sigma <- 2
N <- 1e6

DGP <- function(N, mean=mu, sd=sigma) {
   x <- rnorm(N, mean, sd)
   x
}

loglik <- function(par) {
   mu <- par[1]
   sigma <- par[2]
   sum(dnorm(x, mu, sigma, log=TRUE))
}

x <- DGP(N)
comm.cat("mean x", mean(x), "\n", all.rank=TRUE)

finalize()
