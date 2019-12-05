

#' @export
pNBIIot <- NULL

#' @export
dNBIIot <- NULL

#' @export
qNBIIot <- NULL

#' @export
rNBIIot <- NULL

#' @export
NBIIot <- NULL


#' @export
pNBIIzt <- NULL

#' @export
dNBIIzt <- NULL

#' @export
qNBIIzt <- NULL

#' @export
rNBIIzt <- NULL

#' @export
NBIIzt <-  NULL

#' @export
#' @import rlang
truncate_dist <- function(distr) {
  prefixes <- c("p", "d", "q", "r")
  pdqr_names <- paste0(prefixes, distr)
  pdqr <- parse_exprs(paste0(prefixes, distr))
  names(pdqr) <- prefixes
  e <- caller_env()
  with(pdqr, {
    function_env <- new_environment(data = mget(pdqr_names, envir = e, inherits = TRUE), parent = e)
    dtr <- rlang::new_function(
      args = c(eval_bare(expr(formals(!!d)), env = e), exprs(lower = -Inf, upper = Inf)), 
      body = expr({
        d_args <- formals(sys.function(sys.parent(n = 0)))
        d_args_supplied <- as.list(match.call()[-1])
        d_args[names(d_args_supplied)] <- d_args_supplied
        d_args$lower <- NULL
        d_args$upper <- NULL
        d <- do.call(!!d, d_args)
        
        p_args <- formals(!!p)
        p_args_supplied <- as.list(match.call()[-1])
        p_args_supplied <- p_args_supplied[names(p_args_supplied) %in% names(p_args)]
        p_args[names(p_args_supplied)] <- p_args_supplied
        p_args$q <- NULL    
        p_args$log.p <- FALSE
        
        p_lower <- do.call(!!p, c(list(q = lower), p_args))
        p_upper <- do.call(!!p, c(list(q = upper), p_args))
        p_within <- pmax(p_upper - p_lower, 0)
        d <- if(log) d - log(p_within) else d / p_within
        d <- d * (x>lower) * (x<=upper)
        #if(any(is.nan(d))) browser()
        d
      }), 
      env = function_env
    )
    
    ptr <- rlang::new_function(
      args = c(eval_bare(expr(formals(!!p)), env = e), exprs(lower = -Inf, upper = Inf)), 
      body = expr({
        p_args <- formals(sys.function(sys.parent(n = 0)))
        p_args_supplied <- as.list(match.call()[-1])
        p_args[names(p_args_supplied)] <- p_args_supplied
        p_args$lower <- NULL
        p_args$upper <- NULL
        p_args$q <- NULL    
        p <- do.call(!!p, c(list(q = q), p_args))
        
        p_args$log.p <- FALSE
        p_lower <- do.call(!!p, c(list(q = lower), p_args))
        p_upper <- do.call(!!p, c(list(q = upper), p_args))
        p_within <- pmax(p_upper - p_lower, 0)
        p <- if(log.p) log(exp(p)-p_lower) - log(p_within) else (p-p_lower) / p_within
        p[q<=lower] <- 0
        p[q>upper] <- 1
        p[upper<lower] <- NaN[length(p)>0]
        #if(any(is.nan(p))) browser()
        p
      }), 
      env = function_env
    )
    assign(paste0("d", "t", distr), dtr, envir = e)
    assign(paste0("p", "t", distr), ptr, envir = e)
    invisible(list(dtr = dtr, ptr = ptr))
  })
} 

# ------------------

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

#' prediction function for ggplot2's stat_smooth
#' 
#' @param model Input model
#' @param xseq x-values used for prediction
#' @param se Predict error or not
#' @param level Confidence level
#' @export
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
isoreg_up <-  function(formula, data, ...) {
  M <- model.frame(formula, data)
  f <- as.stepfun(isoreg(M[,2], -M[,1]))
  f2 <- function(x) -f(x)
  class(f2) <- class(f)
  f2
}

#' @export
isoreg_down <-  function(formula, data, ...) {
  M <- model.frame(formula, data)
  as.stepfun(isoreg(M[,2], M[,1]))
}

#' @export
predict.stepfun <- function(model, newdata, se.fit, level, interval) model(newdata[[1]])
