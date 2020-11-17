#' @import gamlss.dist
#' @import gamlss.tr
.onLoad <- function(libname, pkgname) {
  nsname <-"gamlss.dist"
  attname <- paste0("package:", nsname)
  if (!(attname %in% search())) {
    attachNamespace(nsname)
    on.exit(detach(attname, character.only = TRUE))
  }

  fam <- "NBII"
  pars <- list(zt=c(0), ot=c(1), dt=c(1,10))
  types<-list(zt="left", ot="left")
  for (name in c("zt", "ot")) {
    gen_name <- paste0(name, "_")
    # generate truncated negative binomial
    utils::capture.output(eval(expr(gamlss.tr::gen.trun(par=pars[[name]], family=!!fam, name=!!gen_name, type=types[[name]]))))
    for (f in c("p", "d", "r", "q", "")) {
      f_ori_name <- paste0(f, fam)
      f_gen_name <- paste0(f, fam, gen_name)
      f_gen_sym <- parse_expr(f_gen_name)
      f_new_name <- paste0(f, fam, name)

      f_new_formals <- formals(f_ori_name)
      f_args <- parse_exprs(names(f_new_formals))
      names(f_args) <- names(f_new_formals)
      new_f <- local({
        f_gen_name <- f_gen_name
        eval(expr(!!f_gen_sym <- !!f_gen_sym))
        new_f <- function() {}
        formals(new_f) <- f_new_formals
        body(new_f) <- expr({
          args <- formals(sys.function(sys.parent(n = 0)))
          args2 <- as.list(match.call()[-1])
          args[names(args2)] <- args2
          do.call(!!f_gen_sym, args, envir = parent.frame())
        })
        new_f})
      utils::assignInMyNamespace(f_new_name, new_f)
      #assign(f_new_name, new_f, envir= global_env())
      #if(exists(f_gen_name, envir=sys.frame(0), mode="function", inherits=FALSE))
      #  rm(list = f_gen_name, envir = sys.frame(0))
    }
  }
}
