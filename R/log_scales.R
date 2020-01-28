
log_breaks = function(maj, radix=10) {
  function(x) {
    minx         = floor(min(logb(x,radix), na.rm=T)) - 1
    maxx         = ceiling(max(logb(x,radix), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    if (maj) {
      breaks = major_breaks
    } else {
      steps = logb(1:(radix-1),radix)
      breaks = rep(steps, times=n_major) +
        rep(major_breaks, each=radix-1)
    }
    radix^breaks
  }
}

#' Nice log / log1p Scale
#' with radix dependent breaks.
#' @import ggplot2
#' @name log_eng
#' @export
scale_x_log_eng = function(..., radix=10) {
  scale_x_continuous(...,
                     trans=scales::log_trans(radix),
                     breaks=log_breaks(TRUE, radix),
                     minor_breaks=log_breaks(FALSE, radix))
}


#' @rdname log_eng
#' @export
scale_x_log <- scale_x_log_eng


#' @export
#' @rdname log_eng
scale_y_log_eng = function(..., radix=10) {
  scale_y_continuous(...,
                     trans=scales::log_trans(radix),
                     breaks=log_breaks(TRUE, radix),
                     minor_breaks=log_breaks(FALSE, radix))
}

#' @rdname log_eng
#' @export
scale_y_log <- scale_y_log_eng


#' @export
#' @rdname log_eng
log1p_breaks <- function(maj, radix=10) {
  function(x) {
    if ( max(x)>1){
      minx         = -1
      maxx         = ceiling(logb(max(x, na.rm=T), radix)) + 1
      n_major      = maxx - minx + 1
      major_breaks = seq(minx, maxx, by=1)
      if (maj) {
        breaks = major_breaks
      } else {
        steps = logb(1:(radix-1), radix)
        breaks = rep(steps, times=n_major) +
          rep(major_breaks, each=radix-1)
      }
      breaks <- radix^breaks
      if (maj) {
        breaks[1] <- 0
      }
      breaks
    } else {
      if (maj) {
        maxx = ceiling(logb(max(x, na.rm=T), radix))
        seq(0, radix^(maxx+1), radix^(maxx-1))
      } else {
        maxx = ceiling(logb(max(x, na.rm=T), radix))
        seq(0, radix^(maxx+1), radix^(maxx-2))
      }
    }
  }
}

log1p_trans <- function(radix) scales::trans_new("log1p", transform = function(x) log1p(x)/log(radix), inverse = function(x) expm1(x*log(radix)))

#' @export
#' @rdname log_eng
scale_x_log1p = function(..., radix=10) {
  scale_x_continuous(...,
                     trans=log1p_trans(radix),
                     breaks=log1p_breaks(TRUE, radix),
                     minor_breaks=log1p_breaks(FALSE, radix))
}


#' @export
#' @rdname log_eng
scale_y_log1p = function(..., radix=10) {
  scale_y_continuous(...,
                     trans=log1p_trans(radix),
                     breaks=log1p_breaks(TRUE, radix),
                     minor_breaks=log1p_breaks(FALSE, radix))
}
