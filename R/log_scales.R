
#' @export
#' @rdname log_eng
format_selected_breaks <- function(breaks_to_label, format) function(x) ifelse(as.character(x) %in% as.character(breaks_to_label(x[!is.na(x)])), format(x), "")

#' @export
#' @rdname log_eng
format_major_breaks <- format_selected_breaks


#' @rdname log_eng
#' @export
log_breaks = function(maj, radix=10, max_major_breaks_to_show_minor=4) {
  function(x, limits, n) {
    if(!missing(limits)) x <- c(x,limits)
    x <- logb(x[x>=0],radix)
    x <- x[is.finite(x)]
    minx         = floor(min(x)) - 1
    maxx         = ceiling(max(x)) + 1
    n_major      = maxx - minx + 1
    if(!is.finite(minx) |!is.finite(maxx)) return(NA)
    major_breaks = seq(minx, maxx, by=1)
    if (maj | (n_major>max_major_breaks_to_show_minor)) {
      breaks = major_breaks
    } else {
      steps = logb(1:(radix-1),radix)
      breaks = rep(steps, times=n_major) +
        rep(major_breaks, each=radix-1)
    }
    radix^breaks
  }
}

#' @export
#' @rdname log_eng
logX_trans <- function(radix, show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE), format = scales::format_format()) scales::trans_new(
  "logX",
  transform = function(x) log(x)/log(radix),
  inverse = function(x) exp(x*log(radix)),
  breaks = log_breaks(!show_minor_ticks, radix),
  minor_breaks = log_breaks(FALSE, radix),
  format = if(show_minor_ticks) format_selected_breaks(log_breaks(TRUE, radix), format = format) else format, #function(x) {this_major_breaks <- log1p_scaled_breaks(TRUE, radix, scale)(x[!is.na(x)]); ifelse(as.character(x) %in% as.character(this_major_breaks), format(x), "")})
  domain = c(0, Inf)
)

make_log_scale <- function (scale_function) {
  # Substitute will get the unevaluated expression from the calling position
  scale_function <- substitute(scale_function)
  # Create the function template
  func_body <- substitute( {
    SCALE_FUNCTION(..., trans=logX_trans(radix,  show_minor_ticks))
  }, list(SCALE_FUNCTION = scale_function))

  func <- eval(quote(function(..., radix=10, show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE)) {}), envir = parent.frame())
  body(func) <- func_body
  func
}

#' Nice log / log1p / log1pscaled scales
#' with radix dependent breaks.
#' @export
#' @rdname log_eng
scale_x_log <- make_log_scale(scale_x_continuous)

#' @export
#' @rdname log_eng
scale_y_log <- make_log_scale(scale_y_continuous)

#' @export
#' @rdname log_eng
scale_color_log <- make_log_scale(scale_color_continuous)

#' @export
#' @rdname log_eng
scale_fill_log <- make_log_scale(scale_fill_continuous)


#' @rdname log_eng
#' @export
scale_x_log_eng <- scale_x_log

#' @rdname log_eng
#' @export
scale_y_log_eng <- scale_y_log






#' @export
#' @rdname log_eng
log1p_scaled_breaks <- function(maj, radix=10, scale = 1) {
  function(x, limits, n) {
    x <- x * scale
    maxx <- max(x, na.rm = TRUE)
    breaks <- if ( maxx >1){
      minx         = -1
      maxx         = floor(logb(maxx, radix)) + 1
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
    } else if(is.finite(maxx)) {
      if (maj) {
        maxx = ceiling(logb(maxx, radix))
        seq(0, radix^(maxx+1), radix^(maxx-1))
      } else {
        maxx = ceiling(logb(maxx,  radix))
        seq(0, radix^(maxx+1), radix^(maxx-2))
      }
    } else return(numeric(0))
    breaks/scale
  }
}

#' @export
#' @rdname log_eng
log1p_scaled_trans <- function(radix, scale, show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE), format = scales::format_format()) scales::trans_new(
  "log1p_scaled",
  transform = function(x) log1p(x*scale)/log(radix),
  inverse = function(x) expm1(x*log(radix))/scale,
  breaks = log1p_scaled_breaks(!show_minor_ticks, radix, scale),
  minor_breaks = log1p_scaled_breaks(FALSE, radix, scale),
  format = if(show_minor_ticks) format_major_breaks(log1p_scaled_breaks(TRUE, radix, scale), format = format) else format, #function(x) {this_major_breaks <- log1p_scaled_breaks(TRUE, radix, scale)(x[!is.na(x)]); ifelse(as.character(x) %in% as.character(this_major_breaks), format(x), "")})
  domain = c(0, Inf)
)
#' @importFrom ggplot2 expansion
make_log1p_scaled_scale <- function (scale_function) {
  # Substitute will get the unevaluated expression from the calling position
  scale_function <- substitute(scale_function)
  # Create the function template
  func_body <- substitute( {
    SCALE_FUNCTION(..., expand = expand, trans=log1p_scaled_trans(radix, scale, show_minor_ticks))
  }, list(SCALE_FUNCTION = scale_function))

  func <- eval(quote(function(..., radix=10, scale = 1, expand = expansion(mult=c(0, 0.05)), show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE)) {}), envir = parent.frame())
  body(func) <- func_body
  func
}


#' @export
#' @rdname log_eng
scale_x_log1p_scaled <- make_log1p_scaled_scale(scale_x_continuous)

#' @export
#' @rdname log_eng
scale_y_log1p_scaled <- make_log1p_scaled_scale(scale_y_continuous)

#' @export
#' @rdname log_eng
scale_color_log1p_scaled <- make_log1p_scaled_scale(scale_color_continuous)

#' @export
#' @rdname log_eng
scale_fill_log1p_scaled <- make_log1p_scaled_scale(scale_fill_continuous)










#' @export
#' @rdname log_eng
log1p_breaks <- function(maj, radix=10) log1p_scaled_breaks(maj, radix)

#' @export
#' @rdname log_eng
log1p_trans <- function(radix = 10, show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE), format = scales::format_format()) {
  trans <- log1p_scaled_trans(radix, scale = 1, show_minor_ticks = show_minor_ticks, format = format)
  trans$name <- "log1p"
  trans
}

make_log1p_scale <- function (scale_function) {
  # Substitute will get the unevaluated expression from the calling position
  scale_function <- substitute(scale_function)
  # Create the function template
  func_body <- substitute( {
      SCALE_FUNCTION(..., expand = expand, trans=log1p_trans(radix, show_minor_ticks))
  }, list(SCALE_FUNCTION = scale_function))

  func <- eval(quote(function(..., radix=10, expand = expansion(mult=c(0, 0.05)),  show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE)) {}), envir = parent.frame())
  body(func) <- func_body
  func
}


#' @export
#' @rdname log_eng
scale_x_log1p <- make_log1p_scale(scale_x_continuous)

#' @export
#' @rdname log_eng
scale_y_log1p <- make_log1p_scale(scale_y_continuous)

#' @export
#' @rdname log_eng
scale_color_log1p <- make_log1p_scale(scale_color_continuous)

#' @export
#' @rdname log_eng
scale_fill_log1p <- make_log1p_scale(scale_fill_continuous)
#' @export
#' @rdname log_eng
scale_color_log1p <- function(..., radix=10, show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE)) scale_color_continuous(..., trans=log1p_trans(radix=radix, show_minor_ticks = show_minor_ticks))

#' @export
#' @rdname log_eng
scale_fill_log1p <- function(..., radix=10, show_minor_ticks = getOption("shelpr.show_minor_ticks", TRUE)) scale_fill_continuous(..., trans=log1p_trans(radix=radix, show_minor_ticks = show_minor_ticks))


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
