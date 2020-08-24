library(data.table)

#' @export
#' @import data.table
adjust_q_vals <- function(q_vals, ...) 
  UseMethod("adjust_q_vals")

adjust_q_vals.default <- function(q_vals, sets) {
  dt <- data.table(q_val = q_vals, set = sets)
  selections_dt <- dt[, .(min_q_val = min(q_val)), by = .(set)][order(min_q_val)][, n_sets := .I][]
  selections_dt[, effective_min_q_val := pmin(pmax(min_q_val, c(0, cummax(min_q_val*max(n_sets)/n_sets)[-.N])), 1)]
  selections_dt
  dt <- selections_dt[, .(effective_min_q_val, n_sets)][dt[, q_val_ := q_val], .(q_val_adj = q_val * max(n_sets)/n_sets, q_val, set), on = .(effective_min_q_val=q_val_), roll = TRUE, mult = "first"]
  dt[, pmin(q_val_adj, 1)]
}



adjust_q_vals.list <- function(q_vals, ...) {
  dt <- data.table(q_val = unlist(q_vals), set = rep(names(q_vals), sapply(q_vals, length)))
  dt[, q_val_adj := adjust_q_vals(q_vals = q_val, sets = set, ...)]
  lapply(split(dt, by = "set"), `[[`, "q_val_adj")
}


if(FALSE) {

q_vals <- list(a = c(0.01, 0.04, 0.5), b = c(0.011, 0.03, 0.6, 0.7), c = c(0.2))
q_vals
q_vals_adj <- adjust_q_vals(q_vals)
q_vals_adj

mean(replicate(1000, {
  n_h0 <- 10
  n_h1 <- 0
  is_h0 <- c(rep(TRUE, n_h0), rep(TRUE, n_h1))
  p_vals <- c(replicate(n_h0, {t.test(rnorm(10))$p.value}), replicate(n_h1, {t.test(rnorm(10, 1))$p.value}))
  q_vals <- p.adjust(p_vals, method = "BH")
  FDR <- if(any(q_vals<0.1)) mean(is_h0[q_vals<0.1]) else 0
  FDR
}))

mean(replicate(1000, {
  n_h0 <- 10
  n_h1 <- 10
  is_h0 <- c(rep(TRUE, n_h0), rep(FALSE, n_h1))
  p_vals <- c(replicate(n_h0, {t.test(rnorm(10))$p.value}), replicate(n_h1, {t.test(rnorm(10, 1))$p.value}))
  q_vals <- p.adjust(p_vals, method = "BH")
  FDR <- if(any(q_vals<0.1)) mean(is_h0[q_vals<0.1]) else 0
  FDR
}))



mean(replicate(1000, {
  n_h0_0 <- 1000
  n_h1_0 <- 0
  n_h0_1 <- 1000
  n_h1_1 <- 0
  is_h0 <- rep(c(TRUE, FALSE, TRUE, FALSE), c(n_h0_0, n_h1_0,n_h0_1,n_h1_1))
  p_vals_0 <- c(replicate(n_h0_0, {t.test(rnorm(10))$p.value}), replicate(n_h1_0, {t.test(rnorm(10, 1))$p.value}))
  p_vals_1 <- c(replicate(n_h0_1, {t.test(rnorm(10))$p.value}), replicate(n_h1_1, {t.test(rnorm(10, 1))$p.value}))
  q_vals_0 <- p.adjust(p_vals_0, method = "BH")
  q_vals_1 <- p.adjust(p_vals_1, method = "BH")
  q_vals <- c(q_vals_0, q_vals_1)
  FDR <- if(any(q_vals<0.1)) mean(is_h0[q_vals<0.1]) else 0
  FDR
}))

mean(replicate(10000, {
  n_h0_0 <- 10000
  n_h1_0 <- 100
  n_h0_1 <- 10000
  n_h1_1 <- 300
  is_h0 <- rep(c(TRUE, FALSE, TRUE, FALSE), c(n_h0_0, n_h1_0,n_h0_1,n_h1_1))
  p_vals_0 <- c(runif(n_h0_0), runif(n_h1_0)^2)
  p_vals_1 <- c(runif(n_h0_1), runif(n_h1_1)^2)
  q_vals_0 <- p.adjust(p_vals_0, method = "BH")
  q_vals_1 <- p.adjust(p_vals_1, method = "BH")
  q_vals <- c(q_vals_0, q_vals_1)
  FDR <- if(any(q_vals<0.1)) mean(is_h0[q_vals<0.1]) else 0
  FDR
}))

mean(replicate(1000, {
  n_h0_0 <- 10
  n_h1_0 <- 2
  n_h0_1 <- 10
  n_h1_1 <- 0
  is_h0 <- rep(c(TRUE, FALSE, TRUE, FALSE), c(n_h0_0, n_h1_0,n_h0_1,n_h1_1))
  p_vals_0 <- c(replicate(n_h0_0, {t.test(rnorm(10))$p.value}), replicate(n_h1_0, {t.test(rnorm(10, 1))$p.value}))
  p_vals_1 <- c(replicate(n_h0_1, {t.test(rnorm(10))$p.value}), replicate(n_h1_1, {t.test(rnorm(10, 1))$p.value}))
  q_vals_0 <- p.adjust(p_vals_0, method = "BH")
  q_vals_1 <- p.adjust(p_vals_1, method = "BH")
  q_vals <- c(q_vals_0, q_vals_1)
  FDR <- if(any(q_vals<0.1)) mean(is_h0[q_vals<0.1]) else 0
  FDR
}))



mean(replicate(1000, {
  n_h0_0 <- 10
  n_h1_0 <- 0
  n_h0_1 <- 10
  n_h1_1 <- 0
  is_h0 <- rep(c(TRUE, FALSE, TRUE, FALSE), c(n_h0_0, n_h1_0,n_h0_1,n_h1_1))
  p_vals_0 <- c(replicate(n_h0_0, {t.test(rnorm(10))$p.value}), replicate(n_h1_0, {t.test(rnorm(10, 1))$p.value}))
  p_vals_1 <- c(replicate(n_h0_1, {t.test(rnorm(10))$p.value}), replicate(n_h1_1, {t.test(rnorm(10, 1))$p.value}))
  q_vals_0 <- p.adjust(p_vals_0, method = "BH")
  q_vals_1 <- p.adjust(p_vals_1, method = "BH")
  q_vals <- unlist(adjust_q_vals(list(s0 = q_vals_0, s1 = q_vals_1)))
  FDR <- if(any(q_vals<0.1)) mean(is_h0[q_vals<0.1]) else 0
  FDR
}))


mean(replicate(10000, {
  n_h0_0 <- 0
  n_h1_0 <- 10
  n_h0_1 <- 100
  n_h1_1 <- 0
  s <- rep(c(0L, 1L), c(n_h0_0+n_h1_0, n_h0_1+n_h1_1))
  is_h0 <- rep(c(TRUE, FALSE, TRUE, FALSE), c(n_h0_0, n_h1_0,n_h0_1,n_h1_1))
  p_vals_0 <- c(runif(n_h0_0), runif(n_h1_0)^2)
  p_vals_1 <- c(runif(n_h0_1), runif(n_h1_1)^2)
  q_vals_0 <- p.adjust(p_vals_0, method = "BH")
  q_vals_1 <- p.adjust(p_vals_1, method = "BH")
  q_vals <- c(q_vals_0, q_vals_1)
  q_vals <- adjust_q_vals(q_vals, s)
  FDR <- sum(is_h0[q_vals<0.1]) / pmax(1, sum(q_vals<0.1))
  FDR
}))
}
