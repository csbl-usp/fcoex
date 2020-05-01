context("fcoex methods")


#### Setting variables #####
ncells <- 10
my_counts_df_gene_a <- data.frame(matrix(c(rpois(100,15),rpois(300, 5)), ncol = 10))
my_counts_df_gene_b <- data.frame(matrix(c(rpois(120,15),rpois(280, 5)), ncol = 10))

my_counts_df <- rbind(my_counts_df_gene_a, my_counts_df_gene_b)

number_of_rows = nrow(my_counts_df)
rownames(my_counts_df) = paste0(rep("Gene-", number_of_rows), c(1:number_of_rows))
  
  
target <- as.factor(c(rep("A",5), rep("B",5)))


fc <- new("fcoex", expression=my_counts_df, target)
fc <- discretize(fc)
fc <- find_cbf_modules(fc, n_genes_selected_in_first_step = 40)



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
