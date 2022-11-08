library(dlnm)
library(data.table)
library(dlmtree)
library(mgcv)
# source("sim/ts.dlnm.R")
combine_models <- function(mlist) {
  out <- mlist[[1]]
  iter <- out$mcmcIter
  out$DLM <- do.call(rbind, lapply(1:length(mlist), function(i) {
    d <- mlist[[i]]$DLM
    d$Iter <- d$Iter + (i - 1) * iter
    d
  }))
  out$mcmcIter <- mlist[[1]]$mcmcIter * length(mlist)
  out$nIter <- mlist[[1]]$nIter * length(mlist)
  colnames(out$DLM) <- colnames(mlist[[1]]$DLM)
  out$sigma2 <- do.call(c, lapply(mlist, function(l) l$sigma2))
  out$kappa <- do.call(c, lapply(mlist, function(l) l$kappa))
  out$nu <- do.call(c, lapply(mlist, function(l) l$nu))
  out$tau <- do.call(rbind, lapply(mlist, function(l) l$tau))
  out$termNodes <- do.call(rbind, lapply(mlist, function(l) l$termNodes))
  out$gamma <- do.call(rbind, lapply(mlist, function(l) l$gamma))
  out$zirt <- do.call(rbind, lapply(mlist, function(l) l$zirt))
  out$zirtCov <- do.call(rbind, lapply(mlist, function(l) l$zirtCov))
  out$timeProbs <- do.call(rbind, lapply(mlist, function(l) l$timeProbs))
  out$timeCounts <- do.call(rbind, lapply(mlist, function(l) l$timeCounts))
  out$sigma2 <- do.call(c, lapply(mlist, function(l) l$sigma2))
  return(out)
}


##### Sim setup #####
set.seed(sim.num)
nburn <- 2000
niter <- 2000
ntrees <- 20
nthin <- 10
lags <- 20
eval.erc.grid <- 3:30
sim <- list(sim.num = sim.num, n = n, lags = lags, erc = erc, trc = trc,
            mod = list(), summary = list(), metrics <- list())


##### Create data #####
data("chicagoNMMAPS")
setDT(chicagoNMMAPS)
dat <- chicagoNMMAPS[month > 4 & month < 10, .(year, month, dow, temp)]
dat[, paste0("temp.", 0:lags) := shift(temp, 0:lags, type = "lag")]
dat <- dat[sample(which(complete.cases(dat)), n)]
dlnm <- ts.dlnm(erc, trc)
cols <- paste0("temp.", 0:lags)
f <- colSums(sapply(1:n, function(i) dlnm(as.numeric(dat[i, ..cols]), 0:lags)))
truth <- sapply(0:lags, function(t) dlnm(eval.erc.grid, t))
sim$truth <- truth
dat$y <- f + rnorm(n, 0, 1.7 * sd(f))
cor(dat$y, f)
sim[['cor']] <- cor(dat$y, f)
cols <- paste0("temp.", 0:lags)
temp_dat <- as.matrix(dat[, ..cols])
temp_splits <- seq(min(temp_dat), max(temp_dat), length.out = 12)
se <- sd(temp_dat) / 2
cen = 20
out_trans <- function(x) x



##### Monotone-TDLNM #####
m_name <- "tdlnmvs"
mlist <- list()
for (i in 1:restarts) {
  # cat(".")
  set.seed(i)
  mlist[[i]] <- tdlnm(y ~ 1,
                      data = dat,
                      exposure.data = temp_dat,
                      exposure.splits = temp_splits,
                      exposure.se = se,
                      n.trees = ntrees, n.burn = nburn, n.iter = niter, n.thin = nthin,
                      shrinkage = 0,
                      monotone = T, verbose = F)
}
sim$mod[[m_name]] <- combine_models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen, pred.at = eval.erc.grid, verbose = F, conf.level = 0.95)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
sim$metrics[[m_name]] <- list(
  bias = sim$summary[[m_name]]$matfit - truth,
  mse = (sim$summary[[m_name]]$matfit - truth)^2,
  cov = mean((sim$summary[[m_name]]$cilower - 0.01 <= truth) & (sim$summary[[m_name]]$ciupper + 0.01 >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb >= 0.95)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb) >= 0.95) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) > 0))

##### Monotone-TDLNM informative priors #####
m_name <- "tdlnmvs_ip"
truth_lags <- which(colSums(truth > 0) > 0)
zirtp0 <- rep(0.5, lags + 1)
zirtp0[truth_lags] <- 0.9
ts0 <- rep(1, lags)
ts0[c(1, truth_lags)] <- 10

mlist <- list()
for (i in 1:restarts) {
  # cat(".")
  set.seed(i)
  mlist[[i]] <- tdlnm(y ~ 1,
                      data = dat,
                      exposure.data = temp_dat,
                      exposure.splits = temp_splits,
                      exposure.se = se,
                      n.trees = ntrees, n.burn = nburn, n.iter = niter, n.thin = nthin,
                      shrinkage = 0, zirt.p0 = zirtp0, zirt.p0.strength = 2,
                      tree.time.split.params = ts0,
                      monotone = T, verbose = F)
}
sim$mod[[m_name]] <- combine_models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen, pred.at = eval.erc.grid, verbose = F, conf.level = 0.95)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
sim$metrics[[m_name]] <- list(
  bias = sim$summary[[m_name]]$matfit - truth,
  mse = (sim$summary[[m_name]]$matfit - truth)^2,
  cov = mean((sim$summary[[m_name]]$cilower - 0.01 <= truth) & (sim$summary[[m_name]]$ciupper + 0.01 >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb >= 0.95)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb) >= 0.95) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) > 0))


#### TDLNM #####
m_name <- "tdlnm"
mlist <- list()
for (i in 1:restarts) {
  # cat(".")
  set.seed(i)
  mlist[[i]] <- tdlnm(y ~ 1,
                      data = dat,
                      exposure.data = temp_dat,
                      exposure.splits = temp_splits,
                      exposure.se = se,
                      n.trees = ntrees, n.burn = nburn, n.iter = niter, n.thin = nthin,
                      verbose = F)
}
sim$mod[[m_name]] <- combine_models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen, pred.at = eval.erc.grid, verbose = F)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
sim$metrics[[m_name]] <- list(
  bias = (sim$summary[[m_name]]$matfit - truth),
  mse = ((sim$summary[[m_name]]$matfit - truth)^2),
  cov = mean((sim$summary[[m_name]]$cilower <= truth) & (sim$summary[[m_name]]$ciupper >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) > 0))


##### TDLNM informative prior #####
m_name <- "tdlnm_ip"
truth_lags <- which(colSums(truth > 0) > 0)
ts0 <- rep(1, lags)
ts0[c(1,truth_lags)] <- 10
mlist <- list()
for (i in 1:restarts) {
  # cat(".")
  set.seed(i)
  mlist[[i]] <- tdlnm(y ~ 1,
                      data = dat,
                      exposure.data = temp_dat,
                      exposure.splits = temp_splits,
                      exposure.se = se,
                      n.trees = ntrees, n.burn = nburn, n.iter = niter, n.thin = nthin,
                      verbose = F, tree.time.split.params = ts0)
}
sim$mod[[m_name]] <- combine_models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen, pred.at = eval.erc.grid, verbose = F)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
sim$metrics[[m_name]] <- list(
  bias = (sim$summary[[m_name]]$matfit - truth),
  mse = ((sim$summary[[m_name]]$matfit - truth)^2),
  cov = mean((sim$summary[[m_name]]$cilower <= truth) & (sim$summary[[m_name]]$ciupper >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) > 0))



##### GAM #####
m_name <- "gam"
cb <- crossbasis(temp_dat, c(0, lags),
                 arglag = list(fun = "ps", df = 10, intercept = T),
                 argvar = list(fun = "ps", df = 10))
pen <- cbPen(cb)
sim$mod[[m_name]] <- gam(y ~ cb,
                         data = dat, paraPen = list(cb = pen))
sim$summary[[m_name]] <- crosspred(cb, sim$mod[[m_name]], at = eval.erc.grid, bylag = 1, cen = cen)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$matlow <- out_trans(sim$summary[[m_name]]$matlow)
sim$summary[[m_name]]$mathigh <- out_trans(sim$summary[[m_name]]$mathigh)
sim$metrics[[m_name]] <- list(
  bias = (sim$summary[[m_name]]$matfit - truth),
  mse = ((sim$summary[[m_name]]$matfit - truth)^2),
  cov = mean((sim$summary[[m_name]]$matlow <= truth) & (sim$summary[[m_name]]$mathigh >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) > 0)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) == 0)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) > 0)) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) == 0)) / sum(colSums(truth > 0) > 0))


##### GAM informative prior #####
m_name <- "gam_ip"
cb <- crossbasis(temp_dat, c(0, lags),
                 arglag = list(fun = "ps", df = 10, intercept = T),
                 argvar = list(fun = "ps", df = 10))
pen <- cbPen(cb, addSlag = rep(0:1, c(4, 6)))
sim$mod[[m_name]] <- gam(y ~ cb,
                         data = dat, paraPen = list(cb = pen))
sim$summary[[m_name]] <- crosspred(cb, sim$mod[[m_name]], at = eval.erc.grid, bylag = 1, cen = cen)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$matlow <- out_trans(sim$summary[[m_name]]$matlow)
sim$summary[[m_name]]$mathigh <- out_trans(sim$summary[[m_name]]$mathigh)
sim$metrics[[m_name]] <- list(
  bias = (sim$summary[[m_name]]$matfit - truth),
  mse = ((sim$summary[[m_name]]$matfit - truth)^2),
  cov = mean((sim$summary[[m_name]]$matlow <= truth) & (sim$summary[[m_name]]$mathigh >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) > 0)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) == 0)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) > 0)) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$matlow > 0 | sim$summary[[m_name]]$mathigh < 0) == 0)) / sum(colSums(truth > 0) > 0))




sim$mod <- NULL
save(sim, file = paste0("sim.ts.out/", erc, ".", trc, ".", n, ".", sim.num, ".rda"))
