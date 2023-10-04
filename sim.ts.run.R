# install.packages("dlmtree_0.9.0.1.tar.gz",method="source",repo=NULL)
# devtools::install_github("danielmork/dlmtree", ref = "development")
library(parallel)
library(doParallel)
library(foreach)
source("ts.dlnm.R")

nsims <- 100
restarts <- 2
n <- 1000


for (error in c(2, 4, 8)) {
  cl <- makeCluster(48, "PSOCK")
  registerDoParallel(cl)
  getDoParWorkers()

  ret <- foreach(erc = c("exponential", "sublinear", "linear")) %:%
    foreach(trc = c("quadratic", "linear", "piecewise")) %:%
    foreach(sim.num = 1:nsims,
            .errorhandling = "remove",
            .combine = cbind,
            .export = c("restarts", "ts.dlnm", "n", "error"),
            .verbose = TRUE) %dopar% {
              sink("progress.txt", append = TRUE)
              source("sim.ts.R", local = TRUE)
              cat("Complete: error = ", error, "erc =", erc,
                  "trc =", trc, "sim.num =", sim.num, "\n")
            }
  stopCluster(cl)
}
