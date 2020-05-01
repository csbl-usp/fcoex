context("fcoex methods")


#### Setting variables #####
ncells <- 100
my_counts_matrix <- data.frame(matrix(rpois(20000, 5), ncol = ncells))

target <- as.factor(c(rep("A",50), rep("B",50)))
fc <- new("fcoex", expression=my_counts_matrix, target)
fc <- discretize(fc)
fc <- find_cbf_modules(fc, n_genes = 100)


#### Testing ####
test_that("initializer works", {
  expect_is(fc, "fcoex")
})


test_that("discretization works", {
  expect_is(fc@discretized_expression, "data.frame")
})


test_that("module finder works", {
  expect_is(fc@module_list, "list")
})
