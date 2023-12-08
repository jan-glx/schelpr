
#' @export
deming.fit <- function(x, y, noise_ratio = sd(y)/sd(x)) {
  if(missing(noise_ratio) || is.null(noise_ratio)) noise_ratio <- eval(formals(sys.function(0))$noise_ratio) # this is just a complicated way to write `sd(y)/sd(x)`
  delta <-  noise_ratio^2
  x_name <- deparse(substitute(x))

  s_yy <- var(y)
  s_xx <- var(x)
  s_xy <- cov(x, y)
  beta1 <- (s_yy - delta*s_xx + sqrt((s_yy - delta*s_xx)^2 + 4*delta*s_xy^2)) / (2*s_xy)
  beta0 <- mean(y) - beta1 * mean(x)

  res <- c(beta0 = beta0, beta1 = beta1)
  names(res) <- c("(Intercept)", x_name)
  class(res) <- "Deming"
  res
}

#' Deming regression for ggplot stat_smooth and geom_smooth
#' With bootstrapped confidence intervals
#'
#' @param formula formula for the tls
#' @param data Input data.frame
#' @param noise_ratio numeric scalar of the measurement error of LHS over that of RHS if missing or NULL
#' @param ... Unused argument
#' @export
deming <- function(formula, data, R = 100, noise_ratio = NULL, ...){
  ret <- boot::boot(
    data = model.frame(formula, data),
    statistic = function(data, ind) {
      data <- data[ind, ]
      args <- rlang::parse_exprs(colnames(data))
      names(args) <- c("y", "x")
      rlang::eval_tidy(rlang::expr(deming.fit(!!!args, noise_ratio = noise_ratio)), data, env = rlang::current_env())
    },
    R=R
  )
  class(ret) <- c("Deming", class(ret))
  ret
}

# @importFrom ggplot2 predictdf
# @export
#predictdf <- ggplot2:::predictdf

#' prediction function for ggplot2's stat_smooth
#'
#' @param model Input model
#' @param xseq x-values used for prediction
#' @param se Predict error or not
#' @param level Confidence level
# @rawNamespace S3method(ggplot2:::predictdf,Deming)
# @export
predictdf.Deming <- function(model, xseq, se, level) {
  pred <- as.vector(tcrossprod(model$t0, cbind(1, xseq)))
  if(se) {
    preds <- tcrossprod(model$t, cbind(1, xseq))
    data.frame(
      x = xseq,
      y = pred,
      ymin = apply(preds, 2, function(x) quantile(x, probs = (1-level)/2)),
      ymax = apply(preds, 2, function(x) quantile(x, probs = 1-((1-level)/2)))
    )
  } else {
    return(data.frame(x = xseq, y = pred))
  }
}
