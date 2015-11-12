
mu <- 1
sigma <- 2
N <- 1e7

DGP <- function(N, mean=mu, sd=sigma) {
   x <- rnorm(N, mean, sd)
   cat(N, "normals, mean =", mean(x), " std deviation =", sqrt(var(x)), "\n")
   x
}

loglik <- function(par) {
   mu <- par[1]
   sigma <- par[2]
   sum(dnorm(x, mu, sigma, log=TRUE))
}

search1 <- function(mulim=c(-5,5),
                    sigmalim=c(0.001,5),
                    nGrid=100) {
   ## serial computations
   mus <- seq(from=mulim[1], to=mulim[2], length=nGrid)
   sigmas <- seq(from=sigmalim[1], to=sigmalim[2], length=nGrid)
   grid <- as.matrix(expand.grid(mu=mus, sigma=sigmas))
   points <- seq(length=nrow(grid))
   ##
   res <- lapply(points, function(i) loglik(grid[i,]))
   ##
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
   res
}

search2 <- function(mulim=c(-5,5),
                    sigmalim=c(0.001,5),
                    nGrid=100,
                    cores=detectCores()) {
   ## use mclapply (not on windows!)
   library(parallel)
   mus <- seq(from=mulim[1], to=mulim[2], length=nGrid)
   sigmas <- seq(from=sigmalim[1], to=sigmalim[2], length=nGrid)
   grid <- as.matrix(expand.grid(mu=mus, sigma=sigmas))
   points <- seq(length=nrow(grid))
   ##
   res <- mclapply(points, function(i) loglik(grid[i,]), mc.cores=cores)
   ##
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
   res
}

search3 <- function(mulim=c(-5,5),
                    sigmalim=c(0.001,5),
                    nGrid=100) {
   ## use socketCluster
   library(parallel)
   mus <- seq(from=mulim[1], to=mulim[2], length=nGrid)
   sigmas <- seq(from=sigmalim[1], to=sigmalim[2], length=nGrid)
   grid <- as.matrix(expand.grid(mu=mus, sigma=sigmas))
   points <- seq(length=nrow(grid))
   ## -------------------------------------------------
   cl <- makePSOCKcluster(c("localhost", "localhost"),
                          outfile="cluster_log")
   on.exit(stopCluster(cl), add=TRUE)
   print(cl)
   ##
   res <- parLapply(cl, points, function(i) loglik(grid[i,]))
   ## --------------------------------------------------
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
   res
}

search4 <- function(mulim=c(-5,5),
                    sigmalim=c(0.001,5),
                    nGrid=100) {
   ## use socketCluster
   library(parallel)
   mus <- seq(from=mulim[1], to=mulim[2], length=nGrid)
   sigmas <- seq(from=sigmalim[1], to=sigmalim[2], length=nGrid)
   grid <- as.matrix(expand.grid(mu=mus, sigma=sigmas))
   points <- seq(length=nrow(grid))
   ## -------------------------------------------------
   cl <- makePSOCKcluster(c("localhost", "localhost", "localhost", "localhost"),
                          outfile="cluster_log")
   on.exit(stopCluster(cl), add=TRUE)
   print(cl)
   clusterExport(cl, c("x", "loglik"))
   clusterExport(cl, c("grid"), envir=parent.frame(1))
   ##
   res <- clusterApply(cl, points, function(i) loglik(grid[i,]))
   ## --------------------------------------------------
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
   res
}

search5 <- function(mulim=c(-5,5),
                    sigmalim=c(0.001,5),
                    nGrid=100) {
   ## use socketCluster
   library(parallel)
   mus <- seq(from=mulim[1], to=mulim[2], length=nGrid)
   sigmas <- seq(from=sigmalim[1], to=sigmalim[2], length=nGrid)
   grid <- as.matrix(expand.grid(mu=mus, sigma=sigmas))
   points <- seq(length=nrow(grid))
   ## -------------------------------------------------
   eth0 <- system("ifconfig eth0", intern=TRUE)
   eth0 <- grep("inet addr:", eth0, value=TRUE)
   ip <- sub(".*addr:([^ ]+).*", "\\1", eth0)
   cl <- makePSOCKcluster(c(rep("localhost", 4), rep("puyol", 8)),
                          user="otoomet",
                          homogeneous=FALSE,
                          master=ip,
                          outfile="cluster_log.txt")
   on.exit(stopCluster(cl), add=TRUE)
   print(cl)
   clusterExport(cl, c("x", "loglik"))
   ##
   res <- clusterApplyLB(cl, points, function(i) loglik(grid[i,]))
   ## --------------------------------------------------
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
   res
}

search6 <- function(mulim=c(-5,5),
                    sigmalim=c(0.001,5),
                    nGrid=100) {
   ## use socketCluster
   library(parallel)
   library(snow)
   mus <- seq(from=mulim[1], to=mulim[2], length=nGrid)
   sigmas <- seq(from=sigmalim[1], to=sigmalim[2], length=nGrid)
   grid <- as.matrix(expand.grid(mu=mus, sigma=sigmas))
   points <- seq(length=nrow(grid))
   ## -------------------------------------------------
   cl <- makeMPIcluster(4, outfile="cluster_log.txt")
   on.exit(stopCluster(cl), add=TRUE)
   print(length(cl))
   clusterExport(cl, c("x", "loglik"))
   ##
   res <- clusterApplyLB(cl, points, function(i) loglik(grid[i,]))
   ## --------------------------------------------------
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
   res
}
