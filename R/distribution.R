

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
truncate <- function(distr) {
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
