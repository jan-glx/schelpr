
#' @export
isoreg_up <-  function(formula, data, subset, na.action, contrasts = NULL, ..., decreasing = FALSE) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf, contrasts)
  if(decreasing) y <- -y
  ret <- if(ncol(x)==0) isoreg(y) else isoreg(x, y)
  if(decreasing) {
    ret$y <- -ret$y
    ret$yf <- -ret$yf
  }
  as.stepfun(ret)
}

#' @export
isoreg_down <- isoreg_up
formals(isoreg_down)$decreasing <- TRUE

#' @export
predict.stepfun <- function(model, newdata, se.fit, level, interval) model(newdata[[1]])
