ts.dlnm <- function(erc = "exponential", # also, sublinear, exponential
                    trc = "quadratic") # also, piecewise, exponential
{
  library(dlnm)
  library(data.table)
  dat <- chicagoNMMAPS
  setDT(dat)
  dat[, paste0("temp.", 0:20) := shift(temp, 0:20, type = "lag")]

  if (erc == "linear") {
    erc_fun <- function(x) pmax(0, 0.1 * (x - 25))
  } else if (erc == "sublinear") {
    erc_fun <- function(x) ifelse(x > 25, 0.2 * log1p(abs(x - 25)), 0)
  } else if (erc == "exponential") {
    erc_fun <- function(x) 0.2 * pmax(0, exp(0.25 * (x - 25)) - 1)
  } else if (erc == "identity") {
    erc_fun <- function(x) 1
  }


  if (trc == "piecewise") {
    trc_fun <- function(t) ifelse(t < 4, 20, 0)
  } else if (trc == "linear") {
    trc_fun <- function(t) ifelse(t < 6, 6 * (6 - t), 0)
  } else if (trc == "quadratic") {
    trc_fun <- function(t) ifelse(t <= 8, 0.2 * ((t - 8) ^ 2) * (t + 1), 0)
  } else if (trc == "identity") {
    trc_fun <- function(x) 1
  }

  return(function(x, t) as.numeric(erc_fun(x) * trc_fun(t)))
}
