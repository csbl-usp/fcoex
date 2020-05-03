context("fcoex methods")


#### Setting variables #####
ncells <- 10

set.seed("3")
my_counts_df_gene_a <- data.frame(matrix(c(rpois(1000,15),rpois(3000, 5)), ncol = 100))
my_counts_df_gene_b <- data.frame(matrix(c(rpois(1200,15),rpois(2800, 5)), ncol = 100))

my_counts_df <- rbind(my_counts_df_gene_a, my_counts_df_gene_b)

number_of_rows = nrow(my_counts_df)
rownames(my_counts_df) = paste0(rep("Gene-", number_of_rows), c(1:number_of_rows))
  
  
target <- as.factor(c(rep("A",50), rep("B",50)))


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
  
  module_names = c("Gene-7", "Gene-59", "Gene-76")
  
  genes_in_module_for_gene_76 = c("Gene-48", "Gene-50", "Gene-51", "Gene-76", "Gene-74")
  expect_equal(fc@module_list[["Gene-76"]], genes_in_module_for_gene_76)
})
