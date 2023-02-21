
#' @export
scale_color_d3 <- ggsci::scale_color_d3
#' @export
scale_fill_d3 <- ggsci::scale_fill_d3


#' @export
geom_densityjitter <- purrr::partial(ggbeeswarm::geom_quasirandom, method = "pseudorandom")

#' compute confidence intervals for a quantile (e.g. median)
#' See also https://stats.stackexchange.com/questions/122001/confidence-intervals-for-median
#' @export
quantile_cl <- function(y, q=0.5, conf.level = 0.95, na.rm=TRUE) {
  alpha <- 1 - conf.level
  if (na.rm) y <- y[!is.na(y)]
  n <- length(y)
  l <- qbinom(alpha/2, size=n, prob = q)
  u <- 1 + n - l
  ys <- sort.int(c(-Inf, y, Inf), partial = c(1 + l, 1 + u))
  data.frame(
    y = quantile(y, probs = q, na.rm=na.rm, type = 8),
    ymin = ys[1 + l],
    ymax = ys[1 + u]
  )
}

#' @export
median_cl <- function(y, conf.level=0.95, na.rm=TRUE) quantile_cl(y, q=0.5, conf.level=conf.level, na.rm=na.rm)

#' @export
quantile_ci <- function(y, q=0.5, alpha=0.95, na.rm=FALSE) {
  lifecycle::deprecate_soft(when = "0.0.0.9006", what = "quantile_ci()", with = "quantile_cl()")
  quantile_cl(y = y, q = q, conf.level = alpha, na.rm = na.rm)
}


#' @export
fix_plot_limits <- function(p) p + coord_cartesian(xlim=ggplot_build(p)$layout$panel_params[[1]]$x.range, ylim=ggplot_build(p)$layout$panel_params[[1]]$y.range)



#' @export
#' @import ggplot2
#' @importFrom ggforce StatSina
StatSplitSina <- ggproto("StatSplitSina", StatSina, required_aes = c("x", "y", "split"), compute_panel =  function(self, data, ...) {
  new_data <- self$super()$compute_panel(data, ...)
  side_factor <- (-1) ^ (new_data$split)
  new_data$x_diff <- abs(new_data$x_diff) * side_factor
  new_data$xmin <- new_data$x + pmin(0, new_data$width * side_factor)
  new_data$xmax <- new_data$x + pmax(0, new_data$width * side_factor)
  new_data
})

#' @export
#' @import ggplot2
geom_split_sina <- function (mapping = NULL, data = NULL, stat = "splitSina", position = "dodge",
                             ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomPoint,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...))
}

#' @export
#' @import ggplot2
stat_split_sina <- function (mapping = NULL, data = NULL, geom = "Point", position = "dodge",
                             scale = "area", method = "density", bw = "nrd0", kernel = "gaussian",
                             maxwidth = NULL, adjust = 1, bin_limit = 1, binwidth = NULL,
                             bins = NULL, seed = NA, ..., na.rm = FALSE, show.legend = NA,
                             inherit.aes = TRUE)
{
  method <- match.arg(method, c("density", "counts"))
  layer(data = data, mapping = mapping, stat = StatSplitSina, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(scale = scale, method = method, bw = bw,
                      kernel = kernel, maxwidth = maxwidth, adjust = adjust,
                      bin_limit = bin_limit, binwidth = binwidth, bins = bins,
                      seed = seed, na.rm = na.rm, ...))
}



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



#' @export
geom_cell <- function(..., size = 0.2, raster.dpi = 600) ggrastr::rasterize(ggplot2::geom_point(..., size = size), dpi = raster.dpi)

#' @export
print_by <- function(x, by, plot) {
  x <- as.data.table(x)
  split(x, by=by) %>% lapply(function(dt) print(plot(dt)))
  invisible(x)
}

#' modified after ggplot2::calculate_ellipse
#'
#' MIT License
#' Copyright (c) 2020 ggplot2 authors
#' Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#' The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#'
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
calculate_ellipse_weighted <- function(data, vars, type, level, segments){
   dfn <- 2
  dfd <- nrow(data) - 1

  if (!type %in% c("t", "norm", "euclid")) {
    cli::cli_inform("Unrecognized ellipse type")
    ellipse <- matrix(NA_real_, ncol = 2)
  } else if (dfd < 3) {
    cli::cli_inform("Too few points to calculate an ellipse")
    ellipse <- matrix(NA_real_, ncol = 2)
  } else {
    if (type == "t") {
      v <- MASS::cov.trob(data[,vars[1:2]], wt = data[, vars[3]])
    } else if (type == "norm") {
      v <- stats::cov.wt(data[,vars[1:2]], wt = data[,vars[3]])
    } else if (type == "euclid") {
      v <- stats::cov.wt(data[,vars[1:2]], wt = data[, vars[3]])
      v$cov <- diag(rep(min(diag(v$cov)), 2))
    }
    shape <- v$cov
    center <- v$center
    chol_decomp <- chol(shape)
    if (type == "euclid") {
      radius <- level/max(chol_decomp)
    } else {
      radius <- sqrt(dfn * stats::qf(level, dfn, dfd))
    }
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    ellipse <- t(center + radius * t(unit.circle %*% chol_decomp))
  }

  colnames(ellipse) <- vars[1:2]
  ggplot2:::mat_2_df(ellipse)
}

#' @export
StatEllipseWeighted <- ggproto("StatEllipseWeighted", StatEllipse,
                               required_aes = c("x", "y", "weight"),
                               dropped_aes = c("weight"),
                               compute_group = function(data, scales, type = "t", level = 0.95, segments = 51, na.rm = FALSE) {
                                 calculate_ellipse_weighted(data = data, vars = c("x", "y", "weight"), type = type, level = level, segments = segments)
                               }
)

#' @md
#' @seealso [ggplot2::stat_ellipse()]
#' @export
stat_ellipse_weighted <- function (mapping = NULL, data = NULL, geom = "path", position = "identity",
                                   ..., type = "t", level = 0.95, segments = 51, na.rm = FALSE,
                                   show.legend = NA, inherit.aes = TRUE)
{
  layer(data = data, mapping = mapping, stat = StatEllipseWeighted,
        geom = geom, position = position, show.legend = show.legend,
        inherit.aes = inherit.aes, params = rlang::list2(type = type,
                                                         level = level, segments = segments, na.rm = na.rm,
                                                         ...))
}
