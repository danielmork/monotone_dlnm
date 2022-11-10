rm(list=ls())
gc()
library(tidyr)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(tables)
source("ts.dlnm.R")

restarts <- 20
nsims <- 100

cl <- makeCluster(47, "FORK")
registerDoParallel(cl)
getDoParWorkers()

ret <- foreach(erc = c("linear", "exponential", "sublinear")) %:%
  foreach(trc = c("quadratic", "exponential", "piecewise")) %:%
  foreach(n = c(500)) %:%
  foreach(sim.num = 1:nsims, 
          .errorhandling = "remove", 
          .combine = cbind,
          .export = c("restarts", "ts.dlnm"),
          .verbose = FALSE) %dopar% source("sim.ts.R", local = TRUE)
stopCluster(cl)


plot_dat <- data.frame()
res <- data.frame()
for (erc in c("linear", "exponential", "sublinear")) {
  for (trc in c("quadratic", "exponential", "piecewise")) {
    for (n in c(500)) {
      cat(erc, trc, n, "\n")
      cor <- rep(0, nsims)
      mse <- bias <- cov <- tp <- fp <- ci.width <- list()
      for (sim.num in 1:nsims) {
        load(file = paste0("sim.ts.out/", erc, ".", trc, ".", n, ".", sim.num, ".rda"))
        dlnm <- ts.dlnm(erc, trc)
        truth <- sapply(0:20, function(t) dlnm(3:30, t))


        for (mod in names(sim$metrics)) {
          if (sim.num == 1) {
            tp[[mod]] <- fp[[mod]] <- bias[[mod]] <- cov[[mod]] <- mse[[mod]] <- ci.width[[mod]] <- list()
          }
          bias[[mod]][[sim.num]] <- sim$metrics[[mod]]$bias
          mse[[mod]][[sim.num]] <- sim$metrics[[mod]]$mse
          cov[[mod]][[sim.num]] <- sim$metrics[[mod]]$cov
          tp[[mod]][[sim.num]] <- sim$metrics[[mod]]$tp
          fp[[mod]][[sim.num]] <- sim$metrics[[mod]]$fp
          if (mod == "gam" | mod == "gam_ip") {
            ci.width[[mod]][[sim.num]] <- mean(sim$summary[[mod]]$mathigh - sim$summary[[mod]]$matlow)
          } else {
            ci.width[[mod]][[sim.num]] <- mean(sim$summary[[mod]]$ciupper - sim$summary[[mod]]$cilower)
          }
          cor[sim.num] <- sim$cor

          plot_dat <- rbind.data.frame(
            plot_dat,
            data.frame(lag = rep(0:20, each = length(3:30)),
                       x = rep(3:30, 21),
                       fit = as.numeric(sim$summary[[mod]]$matfit),
                       truth = as.numeric(truth),
                       model = mod,
                       n = n,
                       erc = erc,
                       trc = trc,
                       group = sim.num)
          )
        }
      }
      print(summary(cor))

      res <- rbind.data.frame(
        res,
        cbind.data.frame(
          erc = erc, trc = trc, n = n,
          model = names(mse),
          rmse = sapply(names(mse), function(mod) {
            mean(sqrt(Reduce("+", mse[[mod]]) / length(mse[[mod]])))
          }),
          bias = sapply(names(mse), function(mod) {
            mean(Reduce("+", bias[[mod]]) / length(bias[[mod]]))
          }),
          cov = sapply(names(mse), function(mod) {
            mean(Reduce("+", cov[[mod]]) / length(cov[[mod]]))
          }),
          tp = sapply(names(mse), function(mod) {
            mean(Reduce("+", tp[[mod]])) / length(tp[[mod]])
          }),
          fp = sapply(names(mse), function(mod) {
            mean(Reduce("+", fp[[mod]])) / length(fp[[mod]])
          }),
          ci.width = sapply(names(mse), function(mod) {
            mean(Reduce("+", ci.width[[mod]])) / length(ci.width[[mod]])
          })
        )
      )
    }
  }
}
save(res, file = "sim_results.rda")

res <- res %>%
  mutate(prec = tp / (tp + fp),
         f1 = 2 * tp * (tp / (tp + fp)) / (tp + tp / (tp + fp)))
# toLatex(
tabular(Factor(erc) * Factor(trc) ~
          Format(digits = 2) * (rmse + prec) * Factor(model) * mean, data = res)
# )
# toLatex(
tabular(Factor(erc) * Factor(trc) ~
          Format(digits = 2) * (cov + ci.width) * Factor(model) * mean, data = res)
# )
tabular(Factor(erc) * Factor(trc) ~
            Format(digits = 2) * (tp + fp) * Factor(model) * mean, data = res)


plot_x <- 30
plot_dat %>%
  filter(x == plot_x) %>%
  group_by(lag, model, erc, trc) %>%
  summarize(fit = mean(fit), truth = mean(truth)) %>%
  ggplot(aes(x = lag, y = fit, color = erc, linetype = erc)) +
  geom_line(aes(y = truth), color = "grey", size = 2) +
  geom_line(size = 2) +
  facet_grid(model ~ trc) +
  theme_bw(base_size = 16)

plot_lag <- 2
plot_dat %>%
  filter(lag == plot_lag) %>%
  group_by(x, model, erc, trc) %>%
  summarize(fit = mean(fit), truth = mean(truth)) %>%
  ggplot(aes(x = x, y = fit, color = trc, linetype = trc)) +
  geom_line(aes(y = truth), color = "grey", size = 2) +
  geom_line(size = 2) +
  facet_grid(model ~ erc) +
  theme_bw(base_size = 16)
