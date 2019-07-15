#' find_cbf_modules
#'
#' find_cbf_modules uses Symmetrical Uncertainty as a correlation measure and the FCBF algorithm to
#'
#' 1 - Filter the gene list by correlations to a class (Step 1)
#'
#' and
#'
#' 2 - Determine soft thresholds for coexpression to genes predominantly correlated to a class.
#'
#' @param fc A fcoex object containing a discretized expression table
#' @param FCBF_threshold A threshold for the minimum correlation (as determined by symettrical uncertainty)
#' between each variable and the class used for wrapped FCBF function. Defaults to 0.1.
#' @param verbose Adds verbosity. Defaults to FALSE.
#' @param output Form of the output. For an adjacency matrix, "matrix", for the list of modules, "modules".
#' @return Returns a list with the CBF modules found or a adjacency matrix of the graph
#' @import dplyr
#' @import progress
#' @export
#' @examples
#'  data(scDengue)
#'  exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#'  discrete_expression <- as.data.frame(discretize_exprs(exprs))
#'  head(discrete_expression[,1:4])
#'  infection <- SummarizedExperiment::colData(scDengue)
#'  target <- infection$infection
#'  cbf_modules <- find_modules(discrete_expression,target, FCBF_threshold = 0.06)


find_cbf_modules <- function(discretized_exprs, target, FCBF_threshold = 0.1, output ="modules",verbose = FALSE){
  
  # get the SU scores for each gene
  message('Getting SU scores')
  
  su_ic_vector <- FCBF::get_su(discretized_exprs, target)
  su_ic_vector$gene <- rownames(su_ic_vector)
  
  colnames(su_ic_vector)[1] <- 'SU'
  message('Running FCBF to find module headers')
  fcbf_filtered <- FCBF::fcbf(discretized_exprs, target, thresh = FCBF_threshold, verbose = verbose)
  fcbf_filtered$gene <- rownames(fcbf_filtered)
  
  # R does not like points. Subs for -.
  FCBF_genes <- gsub('\\.', '-', fcbf_filtered$gene)
  
  # get only those with an SU score above a threshold
  SU_threshold <- FCBF_threshold
  
  su_ic_vector_small <- su_ic_vector[su_ic_vector$SU > SU_threshold,]
  
  
  SU_genes <- gsub('\\.', '-',su_ic_vector_small[,2])
  exprs_small <- discretized_exprs[SU_genes ,]
  
  
  # get and adjacency matrix for gene to gene correlation
  su_i_j_matrix <- data.frame(genes =  SU_genes)
  message('Calculating adjacency matrix')
  pb_findclusters <- progress_bar$new(total = length(SU_genes),
                                      format =   "[:bar] :percent eta: :eta")
  # this can surely be improved for speed.
  for (i in SU_genes) {
    pb_findclusters$tick()
    gene_i <- as.factor(discretized_exprs[i, ])
    gene_i_correlates <- FCBF::get_su(x = exprs_small, y = as.factor(exprs_small[i, ]))
    gene_i_correlates$gene <- gsub('\\.', '-',rownames(gene_i_correlates))
    gene_i_correlates <- gene_i_correlates[match(su_i_j_matrix$genes,gene_i_correlates$gene),]
    colnames(gene_i_correlates)[1] <- i
    su_i_j_matrix[, i] <- gene_i_correlates[,1]
    
  }
  filtered_su_i_j_matrix <- data.frame(genes =  SU_genes)
  
  message('Getting modules from adjacency matrix')
  for (i in colnames(su_i_j_matrix[,-1])){
    if (all(gsub("\\.", "-",su_ic_vector$gene[seq_len(length(su_i_j_matrix$genes))]) == as.character(su_i_j_matrix$genes))){
      
      tf_vector <- su_i_j_matrix[,i] > su_ic_vector$SU[seq_len(length(su_i_j_matrix$genes))]
      
      filtered_su_i_j_matrix[,i] <- su_i_j_matrix[,i] * tf_vector
    }
  }
  
  
  list_of_fcbf_modules <- list()
  for (seed in FCBF_genes){
    module_members <- as.character(filtered_su_i_j_matrix$genes[filtered_su_i_j_matrix[,seed]>0])
    module_members <- module_members[!is.na(module_members)]
    if(length(module_members) > 1) {
      list_of_fcbf_modules[[seed]] <- module_members
    }
  }
  
  if (output == "matrix"){
    return(filtered_su_i_j_matrix)
  }else if (output == "modules"){
    return(list_of_fcbf_modules)
  }
}


