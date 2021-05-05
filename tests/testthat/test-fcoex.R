context("fcoex methods")

#### Setting variables #####
ncells <- 10

set.seed("3")
my_counts_df_gene_a <- data.frame(matrix(c(rpois(500, 12), rpois(1500, 5)), ncol = 100))
my_counts_df_gene_b <- data.frame(matrix(c(rpois(600, 12), rpois(1400, 5)), ncol = 100))

my_counts_df <- rbind(my_counts_df_gene_a, my_counts_df_gene_b)

number_of_rows <- nrow(my_counts_df)
rownames(my_counts_df) <- paste0(rep("Gene-", number_of_rows), c(1:number_of_rows))


target <- as.factor(c(rep("A", 50), rep("B", 50)))


fc <- new("fcoex", expression = my_counts_df, target)
fc <- discretize(fc)
fc <- find_cbf_modules(fc, n_genes_selected_in_first_step = 15, is_parallel = FALSE)



#### Testing ####
test_that("initializer works", {
  expect_is(fc, "fcoex")
})


test_that("discretization works", {
  expect_is(fc@discretized_expression, "data.frame")
})


test_that("module finder works", {
  expect_is(fc@module_list, "list")

  genes_in_module_for_gene_25 <- c("Gene-25", "Gene-34")
  expect_equal(fc@module_list[["Gene-25"]], genes_in_module_for_gene_25)
})
