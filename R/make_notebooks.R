#'  These ar used by most functions in this file
#' @import rlang
#' @importFrom magrittr `%>%` `%<>%`
NULL

#' @export
add_notebook <- function(object, ...){
  dots <- ensyms(...)
  dependencies <- tidyselect::vars_select(.vars = names(object$notebooks), !!!dots[-1])
  notebook_file <- tidyselect::vars_select(.vars = sapply(dots[1], as.character), !!!dots[1]) #`mode<-`(dots[1], "character")# tidyselect::vars_select(.vars = sapply(dots[1], as.character), !!!dots[1])
  
  object$notebooks %<>% append(notebook_file)
  object$dependencies %<>% append(list2(!!names(notebook_file):=dependencies))
  object
}

#a <- add_notebook

#' @export
render_separately <- function(...)  callr::r(function(...) rmarkdown::render(..., envir = globalenv()), args = list(...), show = TRUE)

#' @export
knit_all_formats <- function(rmd_file, params_=list(), figure_formats=c("png", "pdf", "svg")) {
  rmd_file_plain <- fs::path_ext_remove(fs::path_file(rmd_file))
  
  params_string <- gsub(" +", " ", paste(deparse(params_[seq_along(params_)]), collapse=""))
  params_hash <- substr(digest::digest(params_string), 1, 8)
  results_dir <- fs::path("results", fs::path_sanitize(rmd_file_plain), fs::path_sanitize(params_hash))
  fs::dir_create(results_dir)
  cat(params_string, file = fs::path(results_dir, "params.txt"))
  cat(params_string, file = fs::path(results_dir,  substr(paste0(fs::path_sanitize(params_string), ext="txt"), 1,50)))
  
  lapply(
    figure_formats,
    function(figure_format) render_separately(
      input = rmd_file,
      params = c(params_, results_dir=results_dir),
      output_file = paste0(rmd_file_plain, "_", figure_format),
      output_format = rmarkdown::html_document(dev=figure_format, keep_md=TRUE),
      output_dir = fs::path(results_dir), 
      encoding = 'UTF-8'
    )
  )
  results_dir
}

hash_params <- function(params) {
  params_string <- gsub(" +", " ", paste(deparse(params[sort(names(params))]), collapse=""))
  params_hash <- substr(digest::digest(params_string), 1, 8)
}

#' @export
bind_parameters <- function(analysis, output_dir, figure_format="png", discard_results = FALSE,...){
  all_params <- as.list(enquos(...))
  unused_params <- all_params
  analysis$commands <- list()
  analysis$out_file <- list()
  for(notebook in names(analysis$notebooks)) {
    notebook_file <- analysis$notebooks[[notebook]]
    params <- rmarkdown::yaml_front_matter(fs::path("notebooks", notebook_file))$params
    params$results_dir <- NULL #results_dir 
    
    augmented_params <- list()
    augmented_params[names(analysis$dependencies[[notebook]])] <- all_params[analysis$dependencies[[notebook]]]
    augmented_params[is.null(augmented_params)] <- NULL
    augmented_params <- c(augmented_params, all_params)
    
    # to inform the user that he might have done wrong:
    params_not_supplied <- setdiff(names(params), names(augmented_params))
    unused_params[names(params)] <- NULL
    if (length(params_not_supplied)>0) message(length(params_not_supplied), " parameter(s) not supplied for \"", notebook, "\". Using defaults:\n", 
                                               paste0(params_not_supplied, ": ", params[params_not_supplied],  collapse="\n"), "\n")
    
    this_notebook_params <- params
    this_notebook_params[names(params)] <- augmented_params[names(params)]
    this_notebook_params_hash <- hash_params(this_notebook_params)
    results_dir <- fs::path(output_dir, fs::path_sanitize(fs::path_ext_remove(notebook_file)), fs::path_sanitize(this_notebook_params_hash))
    this_notebook_params$results_dir <- if(discard_results) tempdir() else results_dir
    
    #save results_dir for dependent
    all_params[notebook] <- results_dir
    out_file <- fs::path(results_dir, paste0(fs::path_ext_remove(notebook_file), "_", figure_format, ".html"))
    analysis$out_file[[notebook]] <- out_file
    ###################################################################################
    analysis$commands[[notebook]] <- expr(rmarkdown::render(
      input = !!(fs::path("notebooks", notebook_file)),
      output_format = rmarkdown::html_document(dev=!!figure_format, keep_md=TRUE),
      output_file = !!out_file,
      output_dir = !!results_dir, 
      params=!!this_notebook_params
      ))
    ####################################################################################
  }
  if(length(unused_params)>0) message(length(unused_params), " parameter(s) supplied but not used: ", paste0(names(unused_params),  collapse=", "))
  analysis
}

expr_to_shell <- function(expr) {
  paste0("Rscript -e '", gsub(pattern = "'", replacement = "'''", paste0(deparse(expr), collapse="")) ,"'")
}

#' @export
generate_makefile <- function(analysis, alt_dependencies = analysis, analysis_name = deparse(substitute(analysis))) {
  substitute_out_files <- analysis$out_file
  substitute_out_files[names(alt_dependencies$out_file)] <- alt_dependencies$out_file
  paste0(analysis_name, ":", paste0(" ", analysis$out_file, collapse=""), "\n") %>% 
    c(sapply(names(analysis$notebooks), function(notebook) { 
      paste0(analysis$out_file[notebook], ": ", fs::path("notebooks", analysis$notebooks[[notebook]]) ,paste0(" ", substitute_out_files[analysis$dependencies[[notebook]]], collapse=""), "\n",
             paste0("\t", expr_to_shell(analysis$commands[[notebook]]), collapse="\n"), "\n")
    })) %>%
    paste0(collapse="\n") %>%
    cat(file=paste0(analysis_name, ".Makefile"))
  if(!fs::file_exists("Makefile")) cat("include *.Makefile\n", file="Makefile")
  invisible(NULL)
}

#' @export
bind_and_make_all_formats <- function(analysis, ..., output_dir = "results", figure_formats = c("png", "svg", "pdf")) {
  analysis_name <- paste0(deparse(substitute(anaysis)), collapse = "")
  params <- enquos(...)
  
  cat(paste0(analysis_name, ":", paste0(" ", analysis_name, "_", figure_formats, collapse=" "), "\n"), file = paste0(analysis_name, ".Makefile"))
  
  main_analysis <- eval_tidy(expr(bind_parameters(analysis, output_dir = output_dir, figure_format = figure_formats[1], discard_results = FALSE, !!!params)))
  generate_makefile(main_analysis, analysis_name = paste0(analysis_name, "_", figure_formats[1]))
  lapply(figure_formats[-1], function(figure_format) {
    expr(bind_parameters(analysis, output_dir = output_dir, figure_format = figure_format, discard_results = TRUE, !!!params)) %>% 
      eval_tidy %>% generate_makefile(main_analysis, analysis_name = paste0(analysis_name, "_", figure_format))
  })
  invisible(NULL)
}

dependency_matrix <- function(dependencies) {
  all_notebooks <- names(dependencies)
  dependencies <- setNames(stack(dependencies), c("parent", "child"))
  n <- length(all_notebooks)
  map <- seq(1, n)
  names(map) <- all_notebooks
  A <- Matrix::sparseMatrix(i = as.integer(map[dependencies$child]), j = as.integer(map[dependencies$parent]), dims=c(n, n), dimnames = list(all_notebooks, all_notebooks))
  A
}

simplify_dependencies <- function(dependencies) {
  A <- dependency_matrix(dependencies)
  n <- ncol(A)
  S <- 0
  An <- A
  for (i in seq_len(n-1)) {
    An <- An %*% A
    S <- S | An
  }
  A & !S
}