library(dlnm)
library(data.table)
library(dlmtree)
library(mgcv)
library(coda)


##### Sim setup #####
set.seed(sim.num)
nburn <- 2000
niter <- 5000
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
dat$y <- f + rnorm(n, 0, error * sd(f))
cor(dat$y, f)
sim[['cor']] <- cor(dat$y, f)
cols <- paste0("temp.", 0:lags)
temp_dat <- as.matrix(dat[, ..cols])
temp_splits <- seq(min(temp_dat), max(temp_dat), length.out = 12)
se <- sd(temp_dat) / 2
cen = 20
out_trans <- function(x) x





##### Monotone-TDLNM #####
m_name <- "m-tdlnm"
mlist <- list()
s0 <- diag(lags + 1) * 0.561^2
for (i in 1:restarts) {
  # cat(".")
  set.seed(i)
  mlist[[i]] <- tdlnm(y ~ 1,
                      data = dat,
                      exposure.data = temp_dat,
                      exposure.splits = temp_splits,
                      exposure.se = se,
                      n.trees = ntrees, n.burn = nburn, n.iter = niter, n.thin = nthin,
                      monotone.sigma = s0, altmcmc = 1,
                      monotone = T, verbose = F)
}
sim$mod[[m_name]] <- combine.models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen,
                                 pred.at = eval.erc.grid, verbose = F,
                                 conf.level = 0.95)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
s <- lapply(mlist, summary, pred.at = eval.erc.grid, verbose = F, mcmc = T)
s_mcmc <- do.call(
  mcmc.list,
  lapply(1:restarts, function(i) {
    d <- as.data.table(as.data.frame.table(s[[i]]$dlm_mcmc))
    mcmc(dcast(d, Var3 ~ Var1 + Var2, value.var = "Freq")[, Var3 := NULL][])
  })
)
sim$summary[[m_name]]$gelman_rubin <- gelman.diag(s_mcmc, autoburnin = F, multivariate = F)

sim$metrics[[m_name]] <- list(
  bias = sim$summary[[m_name]]$matfit - truth,
  mse = (sim$summary[[m_name]]$matfit - truth)^2,
  cov = mean((sim$summary[[m_name]]$cilower - 0.05 <= truth) & (sim$summary[[m_name]]$ciupper + 0.05 >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb >= 0.95)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb) >= 0.95) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) > 0),
  convergence_gr = median(sim$summary[[m_name]]$gelman_rubin$psrf[, 1], na.rm = T))



##### Monotone-TDLNM informative priors #####
m_name <- "m-tdlnm_ip"
truth_lags <- which(colSums(truth > 0) > 0)
g0 <- rep(0, lags + 1) # prior prob 0.005-0.995
g0[truth_lags] <- 4.119 # prior prob 0.95-0.995
s0 <- diag(lags + 1) * 2.701^2
for (s in truth_lags)
  s0[s, s] <- 0.599^2
ts0 <- rep(1, lags)
ts0[truth_lags] <- 10
ts0 <- ts0 / sum(ts0)

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
                      monotone.gamma0 = g0, monotone.sigma = s0,
                      time.split.prob = ts0, altmcmc = 1,
                      monotone = T, verbose = F)
}
sim$mod[[m_name]] <- combine.models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen,
                                 pred.at = eval.erc.grid, verbose = F,
                                 conf.level = 0.95)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
s <- lapply(mlist, summary, pred.at = eval.erc.grid, verbose = F, mcmc = T)
s_mcmc <- do.call(
  mcmc.list,
  lapply(1:restarts, function(i) {
    d <- as.data.table(as.data.frame.table(s[[i]]$dlm_mcmc))
    mcmc(dcast(d, Var3 ~ Var1 + Var2, value.var = "Freq")[, Var3 := NULL][])
  })
)
sim$summary[[m_name]]$gelman_rubin <- gelman.diag(s_mcmc, autoburnin = F, multivariate = F)

sim$metrics[[m_name]] <- list(
  bias = sim$summary[[m_name]]$matfit - truth,
  mse = (sim$summary[[m_name]]$matfit - truth)^2,
  cov = mean((sim$summary[[m_name]]$cilower - 0.05 <= truth) & (sim$summary[[m_name]]$ciupper + 0.05 >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb >= 0.95)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (sim$summary[[m_name]]$splitProb) >= 0.95) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (sim$summary[[m_name]]$splitProb < 0.95)) / sum(colSums(truth > 0) > 0),
  convergence_gr = median(sim$summary[[m_name]]$gelman_rubin$psrf[, 1], na.rm = T))




#### TDLNM #####
m_name <- "tdlnm"
mlist <- list()
for (i in 1:restarts) {
  # cat(".")
  set.seed(i * 100)
  mlist[[i]] <- tdlnm(y ~ 1,
                      data = dat,
                      exposure.data = temp_dat,
                      exposure.splits = temp_splits,
                      exposure.se = se,
                      n.trees = 20, n.burn = nburn, n.iter = niter, n.thin = nthin,
                      verbose = F)
}
sim$mod[[m_name]] <- combine.models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen,
                                 pred.at = eval.erc.grid, verbose = F)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
s <- lapply(mlist, summary, pred.at = eval.erc.grid, verbose = F, mcmc = T)
s_mcmc <- do.call(
  mcmc.list,
  lapply(1:restarts, function(i) {
    d <- as.data.table(as.data.frame.table(s[[i]]$dlm_mcmc))
    mcmc(dcast(d, Var3 ~ Var1 + Var2, value.var = "Freq")[, Var3 := NULL][])
  })
)
sim$summary[[m_name]]$gelman_rubin <- gelman.diag(s_mcmc, autoburnin = F, multivariate = F)

sim$metrics[[m_name]] <- list(
  bias = (sim$summary[[m_name]]$matfit - truth),
  mse = ((sim$summary[[m_name]]$matfit - truth)^2),
  cov = mean((sim$summary[[m_name]]$cilower <= truth) & (sim$summary[[m_name]]$ciupper >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) > 0),
  convergence_gr = median(sim$summary[[m_name]]$gelman_rubin$psrf[, 1], na.rm = T))
#
#
# ##### TDLNM informative prior #####
m_name <- "tdlnm_ip"
truth_lags <- which(colSums(truth > 0) > 0)
ts0 <- rep(1, lags)
ts0[truth_lags] <- 10
ts0 <- ts0 / sum(ts0)
mlist <- list()
for (i in 1:restarts) {
  # cat(".")
  set.seed(i * 100)
  mlist[[i]] <- tdlnm(y ~ 1,
                      data = dat,
                      exposure.data = temp_dat,
                      exposure.splits = temp_splits,
                      exposure.se = se,
                      n.trees = 20, n.burn = nburn, n.iter = niter, n.thin = nthin,
                      time.split.prob = ts0,
                      verbose = F)
}
sim$mod[[m_name]] <- combine.models(mlist)
sim$summary[[m_name]] <- summary(sim$mod[[m_name]], cenval = cen,
                                 pred.at = eval.erc.grid, verbose = F)
sim$summary[[m_name]]$matfit <- out_trans(sim$summary[[m_name]]$matfit)
sim$summary[[m_name]]$cilower <- out_trans(sim$summary[[m_name]]$cilower)
sim$summary[[m_name]]$ciupper <- out_trans(sim$summary[[m_name]]$ciupper)
s <- lapply(mlist, summary, pred.at = eval.erc.grid, verbose = F, mcmc = T)
s_mcmc <- do.call(
  mcmc.list,
  lapply(1:restarts, function(i) {
    d <- as.data.table(as.data.frame.table(s[[i]]$dlm_mcmc))
    mcmc(dcast(d, Var3 ~ Var1 + Var2, value.var = "Freq")[, Var3 := NULL][])
  })
)
sim$summary[[m_name]]$gelman_rubin <- gelman.diag(s_mcmc, autoburnin = F, multivariate = F)

sim$metrics[[m_name]] <- list(
  bias = (sim$summary[[m_name]]$matfit - truth),
  mse = ((sim$summary[[m_name]]$matfit - truth)^2),
  cov = mean((sim$summary[[m_name]]$cilower <= truth) & (sim$summary[[m_name]]$ciupper >= truth)),
  tp = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) > 0),
  tn = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) == 0),
  fp = sum((colSums(truth > 0) == 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) > 0)) / sum(colSums(truth > 0) == 0),
  fn = sum((colSums(truth > 0) > 0) & (colSums(sim$summary[[m_name]]$cilower > 0 | sim$summary[[m_name]]$ciupper < 0) == 0)) / sum(colSums(truth > 0) > 0),
  convergence_gr = median(sim$summary[[m_name]]$gelman_rubin$psrf[, 1], na.rm = T))



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
save(sim, file = paste0("sim.ts.out.newMCMC/", erc, ".", trc, ".", n, ".", error,
                        ".", sim.num, ".rda"))
