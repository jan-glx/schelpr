
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
