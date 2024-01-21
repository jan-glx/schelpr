

partition_funs <- list(find_best_partitioning_non_recursive = find_best_partitioning_non_recursive, find_best_partitioning_recursive = find_best_partitioning_recursive)
for(partition_fun_name in names(partition_funs)) {
  test_that(paste0(partition_fun_name, " works"), {
    partition_fun <- partition_funs[[partition_fun_name]]

    expect_equal(partition_fun(numeric(0), 0), list(partitioning = integer(0), H = 0))
    expect_equal(partition_fun(numeric(0), 1), list(partitioning = integer(0), H = 0))
    expect_equal(partition_fun(numeric(0), 2), list(partitioning = integer(0), H = 0))

    expect_equal(partition_fun(1, 0), list(partitioning = integer(0), H = 0))
    expect_equal(partition_fun(1, 1), list(partitioning = 1L, H = 0))
    expect_equal(partition_fun(1, 2), list(partitioning = 1L, H = 0))

    expect_equal(partition_fun(c(0.5, 0.5), 0), list(partitioning = integer(0), H = 0))
    expect_equal(partition_fun(c(0.5, 0.5), 1), list(partitioning = 2L, H = 0))
    expect_equal(partition_fun(c(0.5, 0.5), 2), list(partitioning = c(1L,2L), H = 1))
    expect_equal(partition_fun(c(0.5, 0.5), 3), list(partitioning = c(1L,2L), H = 1))


    expect_equal(partition_fun(c(0.3, 0.3, 0.4), 0), list(partitioning = integer(0), H = 0))
    expect_equal(partition_fun(c(0.3, 0.3, 0.4), 1), list(partitioning = 3L, H = 0))
    expect_equal(partition_fun(c(0.3, 0.3, 0.4), 2), list(partitioning = c(2L, 3L), H = -(0.4*log2(0.4)+0.6*log2(0.6))))
    expect_equal(partition_fun(c(0.4, 0.3, 0.3), 2), list(partitioning = c(1L, 3L), H = -(0.4*log2(0.4)+0.6*log2(0.6))))
    })
}


test_that(" partionioning function give same results ", {
  withr::with_seed(1, {
    replicate(10, {
      x1 <- sample(1:10,50, prob=1:10, replace=TRUE)
      p1 <- diff(c(0,environment(ecdf(x1))$y))
      p1
      expect_equal(find_best_partitioning_non_recursive(p1, 4), find_best_partitioning_recursive(p1, 4))
    })
  })
})

