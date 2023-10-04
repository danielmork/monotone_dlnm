library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
round2 <- function(...) lapply(list(...), round, digits = 2)


nsims <- 100
n <- 1000

res <- data.table()
for (error in c(2, 4, 8)) {
  for (erc in c("linear", "exponential", "sublinear")) {
    for (trc in c("quadratic", "linear", "piecewise")) {
      cat(erc, trc, n, error, "\n")
      cor <- rep(0, nsims)
      mse <- bias <- cov <- tp <- fp <- ci.width <- gr <- list()
      i <- 1

      for (sim.num in 1:nsims) {
        f <- paste0(erc, ".", trc, ".", n, ".", error,
                    ".", sim.num, ".rda")
        if (!(f %in% list.files("sim.ts.out.newMCMC")))
          next;
        load(file = paste0("sim.ts.out.newMCMC/", f))

        for (mod in names(sim$metrics)) {
          if (i == 1) {
            tp[[mod]] <- fp[[mod]] <- bias[[mod]] <- cov[[mod]] <- list()
            gr[[mod]] <- mse[[mod]] <- ci.width[[mod]] <- list()
          }
          bias[[mod]][[i]] <- sim$metrics[[mod]]$bias
          mse[[mod]][[i]] <- sim$metrics[[mod]]$mse
          cov[[mod]][[i]] <- sim$metrics[[mod]]$cov
          tp[[mod]][[i]] <- sim$metrics[[mod]]$tp
          fp[[mod]][[i]] <- sim$metrics[[mod]]$fp
          if (mod == "gam" | mod == "gam_ip") {
            ci.width[[mod]][[i]] <- mean(sim$summary[[mod]]$mathigh - sim$summary[[mod]]$matlow)
            gr[[mod]][[i]] <- 0
          } else {
            ci.width[[mod]][[i]] <- mean(sim$summary[[mod]]$ciupper - sim$summary[[mod]]$cilower)
            gr[[mod]][[i]] <- sim$metrics[[mod]]$convergence_gr
          }
          cor[i] <- sim$cor
        }
        i <- i + 1
      } # end for 1:nsims
      # print(summary(cor))

      if (length(mse) > 0) {
        res <- rbind(
          res,
          cbind.data.frame(
            erc = erc, trc = trc, n = n, error = error,
            nsims = length(mse[[1]]),
            model = names(mse),
            rmse = sapply(names(mse), function(mod) {
              sqrt(mean(Reduce("+", mse[[mod]]) / length(mse[[mod]]), na.rm = TRUE))
            }),
            bias = sapply(names(mse), function(mod) {
              mean(Reduce("+", bias[[mod]]) / length(bias[[mod]]))
            }),
            bias = sapply(names(bias), function(mod) {
              mean(sapply(1:length(bias), function(n) bias[[mod]][[n]]), na.rm = TRUE)
            }),
            cov = sapply(names(cov), function(mod) {
              mean(sapply(1:length(cov), function(n) cov[[mod]][[n]]), na.rm = TRUE)
            }),
            tp = sapply(names(tp), function(mod) {
              mean(sapply(1:length(tp), function(n) tp[[mod]][[n]]), na.rm = TRUE)
            }),
            fp = sapply(names(fp), function(mod) {
              mean(sapply(1:length(fp), function(n) fp[[mod]][[n]]), na.rm = TRUE)
            }),
            ci.width = sapply(names(ci.width), function(mod) {
              mean(sapply(1:length(ci.width), function(n) ci.width[[mod]][[n]]), na.rm = TRUE)
            }),
            gr = sapply(names(gr), function(mod) {
              mean(sapply(1:length(gr), function(n) gr[[mod]][[n]]), na.rm = TRUE)
            })
          ) # cbind
        ) # rbind
      } # end if length > 1
    } # end for trc
  } # end for erc
} # end for error
res$model2 <- factor(res$model,
                     levels = c("gam", "tdlnm", "m-tdlnm", "gam_ip", "tdlnm_ip", "m-tdlnm_ip"),
                     labels = c("GAM", "TDLNM", "Monotone", "GAM IP", "TDLNM IP", "Monotone IP"))
res$erc <- factor(res$erc, levels = c("exponential", "linear", "sublinear"),
                  labels = c("Exp", "Lin", "Sub"))
res$trc <- factor(res$trc, levels = c("linear", "piecewise", "quadratic"),
                  labels = c("Lin", "Pw", "Quad"))
res[, prec := tp / (tp + fp + 1e-6)]
res[, f1 := 2 * tp * (tp / (tp + fp)) / (tp + tp / (tp + fp))]

dcast(res, error + erc + trc ~ model2, value.var = "nsims", fun = round2)
dcast(res, error + erc + trc ~ model2, value.var = "prec", fun = round2)
dcast(res, error + erc + trc ~ model2, value.var = "tp", fun = round2)
dcast(res, error + erc + trc ~ model2, value.var = "fp", fun = round2)
dcast(res, error + erc + trc ~ model2, value.var = "rmse", fun = round2)
dcast(res, error + erc + trc ~ model2, value.var = "cov", fun = round2)
dcast(res, error + erc + trc ~ model2, value.var = "ci.width", fun = round2)
dcast(res[model != "gam" & model != "gam_ip"],
      error + erc + trc ~ model2, value.var = "gr", fun = round2)
save(res, file = "all_res_newMCMC.rda")
