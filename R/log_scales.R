
#' @rdname log_eng
#' @export
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
scale_x_log_eng  <- function(..., breaks = log_breaks(TRUE, radix), minor_breaks = log_breaks(FALSE, radix), radix=10) {
  scale_x_continuous(...,
                     trans=scales::log_trans(radix),
                     breaks = breaks,
                     minor_breaks = minor_breaks)
}


#' @rdname log_eng
#' @export
scale_x_log <- scale_x_log_eng


#' @export
#' @rdname log_eng
scale_y_log_eng <- function(..., breaks = log_breaks(TRUE, radix), minor_breaks = log_breaks(FALSE, radix), radix=10) {
  scale_y_continuous(...,
                     trans=scales::log_trans(radix),
                     breaks = breaks,
                     minor_breaks = minor_breaks)
}

#' @rdname log_eng
#' @export
scale_y_log <- scale_y_log_eng

#' @export
#' @rdname log_eng
log1p_scaled_breaks <- function(maj, radix=10, scale = 1) {
  function(x, limits, n) {
    x <- x * scale
    breaks <- if ( max(x, na.rm = TRUE)>1){
      minx         = -1
      maxx         = floor(logb(max(x, na.rm=T), radix)) + 1
      n_major      = maxx - minx + 1
      major_breaks = seq(minx , maxx, by=1)
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
      } else {
        breaks <- c(0, breaks)
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
    breaks/scale
  }
}

#' @export
#' @rdname log_eng
format_major_breaks <- function(breaks_to_label, format) function(x) ifelse(as.character(x) %in% as.character(breaks_to_label(x[!is.na(x)])), format(x), "")

#' @export
#' @rdname log_eng
log1p_scaled_trans <- function(radix, scale, format = scales::format_format()) scales::trans_new(
  "log1p_scaled",
  transform = function(x) log1p(x*scale)/log(radix),
  inverse = function(x) expm1(x*log(radix))/scale,
  breaks = log1p_scaled_breaks(FALSE, radix, scale),
  format = format_major_breaks(log1p_scaled_breaks(TRUE, radix, scale), format = format) #function(x) {this_major_breaks <- log1p_scaled_breaks(TRUE, radix, scale)(x[!is.na(x)]); ifelse(as.character(x) %in% as.character(this_major_breaks), format(x), "")})
)

#' @export
#' @rdname log_eng
scale_x_log1p_scaled <- function(..., radix=10, scale = 1) scale_x_continuous(..., trans=log1p_scaled_trans(radix, scale))

#' @export
#' @rdname log_eng
scale_y_log1p_scaled <- function(..., radix=10, scale = 1) scale_y_continuous(..., trans=log1p_scaled_trans(radix, scale))

#' @export
#' @rdname log_eng
scale_color_log1p_scaled <- function(..., radix=10, scale = 1) scale_color_continuous(..., trans=log1p_scaled_trans(radix, scale))




#' @export
#' @rdname log_eng
log1p_breaks <- function(maj, radix=10) log1p_scaled_breaks(maj, radix)

#' @export
#' @rdname log_eng
log1p_trans <- function(radix, format = scales::format_format()) {
  trans <- log1p_scaled_trans(radix, scale = 1, format)
  trans$name <- "log1p"
  trans
}

#' @export
#' @rdname log_eng
scale_x_log1p <- function(..., radix=10) scale_x_continuous(..., trans=log1p_trans(radix))

#' @export
#' @rdname log_eng
scale_y_log1p <- function(..., radix=10) scale_y_continuous(..., trans=log1p_trans(radix))

#' @export
make_quantile_trans <- function(x, format = scales::label_number()) {
  name <- paste0("quantiles_of_", deparse1(substitute(x)))
  xs <- sort(x)
  N <- length(xs)
  transform <- function(x) findInterval(x, xs)/N # find the last element that is smaller
  inverse <- function(q) xs[1+floor(q*(N-1))]

  scales::trans_new(
    name = name,
    transform = transform,
    inverse = inverse,
    breaks =  function(x, n = 5) inverse(scales::extended_breaks()(transform(x), n)),
    minor_breaks = function(x, n = 5) inverse(scales::regular_minor_breaks()(transform(x), n)),
    format = format,
    domain = xs[c(1, N)]
  )
}
