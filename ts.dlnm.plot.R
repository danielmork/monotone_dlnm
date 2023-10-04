source("ts.dlnm.R")
library(ggplot2)

lags <- 20

df <- data.frame()
for (erc in c("linear", "exponential", "sublinear")) {
  dlnm <- ts.dlnm(erc, "identity")
  truth <- sapply(0:lags, function(t) dlnm(15:30, t))
  df <- rbind.data.frame(df, data.frame(x = 15:30,
                                        y = truth[,4],
                                        curve = erc))
}

ggplot(df, aes(x = x, y = y, color = curve, linetype = curve)) +
  geom_line(size = 2) +
  theme_bw(base_size = 32) +
  theme(legend.position = "none", legend.key.width = unit(1, "cm")) +
  labs(x = "x", y = bquote(f[x]), color = "", linetype = "")


df <- data.frame()
for (trc in c("quadratic", "linear", "piecewise")) {
  dlnm <- ts.dlnm("identity", trc)
  truth <- sapply(0:(2*lags)/2, function(t) dlnm(30, t))
  df <- rbind.data.frame(df, data.frame(x = 0:(2*lags)/2,
                                        y = truth,
                                        curve = trc))
}

ggplot(df, aes(x = x, y = y, color = curve, linetype = curve)) +
  geom_line(size = 2) +
  theme_bw(base_size = 32) +
  theme(legend.position = "none", legend.key.width = unit(1, "cm")) +
  labs(x = "Lag", y = bquote(f[l]), color = "", linetype = "")

