
library(pbdMPI, quiet=TRUE)
init()
rank <- comm.rank()
master <- comm.rank() == 0
comm.cat("Job start", as.character(date()), "\n", quiet=TRUE)
##

mu <- 1
sigma <- 2
N <- 1e5
nGrid <- 100
mulim <- c(-5, 5)
sigmalim <- c(0.001, 5)

DGP <- function(N, mean=mu, sd=sigma) {
   x <- rnorm(N, mean, sd)
   x
}

loglik <- function(par) {
   mu <- par[1]
   sigma <- par[2]
   sum(dnorm(x, mu, sigma, log=TRUE))
}

## Now let just the master compute x and broadcast to all workers
x <- 0
                           # need to initialize x for all workers
if(master) {
   x <- DGP(N)
                           # initialize "real" x
}
x <- bcast(x)

## Now lets calculate the grid
mus <- seq(from=mulim[1], to=mulim[2], length=nGrid)
sigmas <- seq(from=sigmalim[1], to=sigmalim[2], length=nGrid)
grid <- as.matrix(expand.grid(mu=mus, sigma=sigmas))
                           # note: every workers has calculated grid
res <- task.pull(seq(length=nrow(grid)), function(i) loglik(grid[i,]))
if(master) {
   res <- matrix(unlist(res), nGrid, nGrid)
   dimnames(res) <- list(mu=formatC(mus, format="f", digits=4, width=7),
                      sigma=formatC(sigmas, format="f", digits=4, width=7))
   ##
   best <- which(res==max(res), arr.ind=TRUE)
   i <- best[1]
   j <- best[2]
   bestMu <- mus[i]
   bestSigma <- sigmas[j]
   cat("Maximum", res[i,j], "at mu =", bestMu, "and sigma =", bestSigma, "\n")
}

finalize()
