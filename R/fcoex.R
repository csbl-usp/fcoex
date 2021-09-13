# Many parts of this code were directly adapted from the
# CEMiTool package. This include chunks of copied and pasted code
# located inside the source code for CEMiTool functions.
# Functions that contained adapted code were explicit denoted.

# Imports ---------

#' @importFrom grDevices rainbow dev.off pdf
#' @importFrom utils write.table head
#' @importFrom stats cutree dist hclust
#' @importFrom methods new 'slot<-' show
#' @importFrom pathwayPCA read_gmt
#' @import SingleCellExperiment
#' @import dplyr
#' @import Matrix
NULL


# Discretize expression ---------

setOldClass("fcoex")

#' Set the discretized expression attribute
#' Uses the discretize_exprs function of the FCBF package
#'
#' @param fc Object of class \code{fcoex}
#' @param method Method applied to all genes for discretization.
#' Methods available: "varying_width"
#'  (Binarization modulated by the number_of_bins param),
#' "mean" (Split in ON/OFF by each gene mean expression),
#' "median" (Split in ON/OFF by each gene median expression),
#' "min_max_\%" (Similat to the "varying width", a binarization threshold
#' in a \% of the min-max range is set. (minmax\% param)),
#' @param number_of_bins Number of equal-width bins for discretization.
#' Note: it is a binary discretization, with the
#' first bin becoming one class ('low') and the other bins, another class
#' ('high').
#' Defaults to 4.
#' @param min_max_cutoff <- Modulator for the "min_max_\%" method.
#' Defaults to 0.25.
#' @return A data frame with the discretized features in the same
#' order as previously
#'
#' @examples
#' library(SingleCellExperiment)
#' data("mini_pbmc3k")
#' targets <- colData(mini_pbmc3k)$clusters
#' exprs <- as.data.frame(assay(mini_pbmc3k, "logcounts"))
#' fc <- new_fcoex(exprs, targets)
#' fc <- discretize(fc)
#' @import FCBF
#' @rdname discretize
#' @export
setGeneric("discretize", function(fc,
                                  number_of_bins = 4,
                                  method = "varying_width",
                                  min_max_cutoff = 0.25) {
  standardGeneric("discretize")
})

#' @rdname discretize
setMethod("discretize", signature("fcoex"),
          function(fc,
                   number_of_bins = 4,
                   method = "varying_width",
                   min_max_cutoff = 0.25) {
            expression_table <- fc@expression
            discretized_expression <-
              FCBF::discretize_exprs(
                expression_table = expression_table,
                number_of_bins = number_of_bins,
                method = method,
                min_max_cutoff = min_max_cutoff
              )
            colnames(discretized_expression) <-
              colnames(expression_table)
            discretized_expression = as.data.frame(ifelse(discretized_expression ==
                                                            "high", 1, 0))
            
            discretized_matrix = as(as.matrix(discretized_expression), "dgCMatrix")
            
            fc@discretized_expression <- discretized_matrix
            
            return(fc)
          })

# Find CBF modules ------

#' find_cbf_modules
#'
#' find_cbf_modules uses Symmetrical Uncertainty as a correlation measure
#'  and the FCBF algorithm to
#' 1 - Filter the gene list by correlations to a class (Step 1)
#' and
#' 2 - Determine soft thresholds for coexpression to genes predominantly
#' correlated to a class.
#'
#' @param fc A fcoex object containing a discretized expression table
#' @param FCBF_threshold A threshold for the minimum correlation (as
#' determined by symettrical uncertainty)
#' between each variable and the class used for wrapped FCBF function.
#' Defaults to 0.1.
#' @param is_parallel Uses package parallel to paralleliza calculations.
#' Defaults to FALSE.
#' @param verbose Adds verbosity. Defaults to TRUE
#' @param n_genes_selected_in_first_step Sets the number of genes to be selected in the first
#' part of the algorithm. If left unchanged, it defaults to NULL and the
#' minimum_su parameter is used.
#' Caution: it overrides the minimum_su parameter altogether.
#'
#' @examples
#' library(SingleCellExperiment)
#' data("mini_pbmc3k")
#' targets <- colData(mini_pbmc3k)$clusters
#' exprs <- as.data.frame(assay(mini_pbmc3k, "logcounts"))
#' fc <- new_fcoex(exprs, targets)
#' fc <- discretize(fc)
#' fc <- find_cbf_modules(fc)
#' @return Returns a list with the CBF modules found or a adjacency matrix of the graph
#' @import dplyr
#' @import parallel
#' @import progress
#' @import FCBF
#' @export
#' @rdname find_cbf_modules
setGeneric("find_cbf_modules", function(fc,
                                        n_genes_selected_in_first_step = NULL,
                                        FCBF_threshold = 0.1,
                                        verbose = TRUE,
                                        is_parallel = FALSE) {
  standardGeneric("find_cbf_modules")
})

#' @rdname find_cbf_modules
setMethod("find_cbf_modules", signature("fcoex"),
          function(fc,
                   n_genes_selected_in_first_step = NULL,
                   FCBF_threshold = 0.1,
                   verbose = TRUE,
                   is_parallel = FALSE) {
            discretized_exprs = as.data.frame(as.matrix(fc@discretized_expression))
            index <- 1:ncol(discretized_exprs)
            discretized_exprs[, index] <-
              lapply(discretized_exprs[, index], factor, levels = c(1, 0))
            
            target <- fc@target
            
            check_rownames(discretized_exprs)
            
            message("Getting SU scores for each gene")
            
            su_to_class <-
              get_su_ranking_in_relation_to_class(discretized_exprs, target)
            
            message("Running FCBF to find module headers")
            
            minimum_su_for_the_fcbf_algorithm <-
              get_minimum_su(n_genes_selected_in_first_step,
                             su_to_class,
                             FCBF_threshold)
            
            module_headers <- run_fcbf_for_module_headers(
              discretized_exprs,
              target,
              n_genes_selected_in_first_step,
              minimum_su = minimum_su_for_the_fcbf_algorithm,
              verbose = FALSE
            )
            
            su_to_class_higher_than_minimum_su <-
              su_to_class[su_to_class[1] > minimum_su_for_the_fcbf_algorithm,]
            
            genes_from_su_ranking <-
              change_dots_for_dashes(su_to_class_higher_than_minimum_su[, 2])
            
            expression_table_only_with_genes_with_high_su <-
              discretized_exprs[genes_from_su_ranking,]
            
            message("Calculating adjacency matrix")
            
            if (!is_parallel) {
              gene_by_gene_su_correlation <-
                get_gene_by_gene_correlation_matrix_in_series(genes_from_su_ranking,
                                                              expression_table_only_with_genes_with_high_su)
            } else {
              gene_by_gene_su_correlation <-
                get_gene_by_gene_correlation_matrix_in_parallel(
                  genes_from_su_ranking,
                  expression_table_only_with_genes_with_high_su,
                  discretized_exprs,
                  su_to_class_higher_than_minimum_su
                )
            }
            
            message("Trimming and getting modules from adjacency matrix")
            
            filtered_gene_by_gene_su_correlation <-
              trim_correlation_matrix(
                genes_from_su_ranking,
                gene_by_gene_su_correlation,
                su_to_class,
                su_to_class_higher_than_minimum_su
              )
            
            list_of_fcbf_modules <-
              get_list_of_modules(module_headers, filtered_gene_by_gene_su_correlation)
            
            message(paste0(as.character(length(list_of_fcbf_modules)), " modules were found."))
            
            fc@selected_genes <- genes_from_su_ranking
            fc@adjacency <- gene_by_gene_su_correlation
            fc@adjacency_trimmed <-
              filtered_gene_by_gene_su_correlation
            fc@module_list <- list_of_fcbf_modules
            return(fc)
          })


## Functions extracted from "find_cbf_modules" ------

#' .get_correlates
#'
#' auxiliary function for find_cbf_modules
#' @param i A gene to be correlated
#' @param gene_by_gene_su_correlation the dataframe with the correlations to be updated
#' @param discretized_exprs the dataframe with discretized expression
#' to extract a gene
#' @param expression_table_only_with_genes_with_high_su the dataframe to after the filtering step
#' @return the updated column of the gene_by_gene_su_correlation
.get_correlates <- function(i,
                            gene_by_gene_su_correlation,
                            discretized_exprs,
                            expression_table_only_with_genes_with_high_su) {
  gene_i <- as.factor(discretized_exprs[i,])
  gene_i_correlates <-
    FCBF::get_su_for_feature_table_and_vector(
      feature_table = expression_table_only_with_genes_with_high_su,
      target_vector = as.factor(expression_table_only_with_genes_with_high_su[i,])
    )
  gene_i_correlates$gene <-
    change_dots_for_dashes(rownames(gene_i_correlates))
  gene_i_correlates <-
    gene_i_correlates[match(gene_by_gene_su_correlation$genes, gene_i_correlates$gene),]
  colnames(gene_i_correlates)[1] <- i
  gene_by_gene_su_correlation[, i] <- gene_i_correlates[, 1]
  gene_by_gene_su_correlation[, i]
}

check_rownames <- function(discretized_exprs) {
  first_name_in_rows <- rownames(discretized_exprs)[1]
  if (first_name_in_rows == "1") {
    stop(
      "The discretized dataframe does not have rownames. That makes fcoex sad! Please, rerun the discretize(fc) method for a expression table with rownames."
    )
  }
  if (any(grepl(" ", rownames(discretized_exprs)))) {
    stop(
      "Oops, there are spaces in at least one of your rownames. That makes fcoex sad! Please, rerun the discretize(fc) method for a expression table with rownames without space."
    )
  }
}

get_su_ranking_in_relation_to_class <-
  function(discretized_exprs, target) {
    su_to_class <-
      FCBF::get_su_for_feature_table_and_vector(discretized_exprs, target)
    su_to_class$gene <-
      change_dots_for_dashes(rownames(su_to_class))
    colnames(su_to_class)[1] <- "SU"
    return(su_to_class)
  }

get_minimum_su <-
  function(n_genes_selected_in_first_step,
           su_to_class,
           minimum_su) {
    # Heuristic to get a minimum_su when user inputs number of genes
    
    if (length(n_genes_selected_in_first_step)) {
      minimum_su <- su_to_class$SU[n_genes_selected_in_first_step]
      return(minimum_su)
    } else {
      return(minimum_su)
    }
  }

# Acessing the fcoex object -------

#' Get the number of modules in a fcoex object
#'
#' This function was copied and adapted from the CEMiTool package.
#'
#' @param fc Object of class \code{fcoex}
#'
#' @return number of modules
#'
#' @rdname nmodules
#' @examples
#' data("fc")
#' nmodules(fc)
#' @export
setGeneric("nmodules", function(fc) {
  standardGeneric("nmodules")
})

#' @rdname nmodules
setMethod("nmodules", signature("fcoex"),
          function(fc) {
            n <- 0
            if ((length(fc@module_list)) > 0) {
              n <- length(fc@module_list)
            } else {
              warning("Run find_cbf_modules function to get modules!")
            }
            return(n)
          })

#' Get the number of genes in modules in a fcoex object
#' This function was copied and adapted from the CEMiTool package.
#'
#' @param fc Object of class \code{fcoex}
#' @param module Default is NULL. If a character string designating a
#' module is#' given, the number of genes in that module is returned instead.
#' @examples
#' data("fc")
#' mod_gene_num(fc, module = "TYROBP")
#' @return The number of genes in module(s)
#'
#' @rdname mod_gene_num
#' @export
setGeneric("mod_gene_num", function(fc, module = NULL) {
  standardGeneric("mod_gene_num")
})
#' @rdname mod_gene_num
setMethod("mod_gene_num", signature(fc = "fcoex"),
          function(fc, module = NULL) {
            if (!is.null(module)) {
              if (!(all(module %in% mod_names(fc)))) {
                stop("Module '", module, "' not in fcoex object!")
              }
            }
            if (!length(module_genes(fc)) > 0) {
              stop("No modules in fcoex object!")
            }
            if (!is.null(module)) {
              mod_genes <- fc@module_list[[module]]
              return(mod_genes)
            }
            
            if (is.null(module)) {
              stop("No modules selected!")
            }
          })

#' Get module names in a fcoex object
#'
#' This function was copied and adapted from the CEMiTool package.
#'
#' @param fc Object of class \code{fcoex}
#' @param include_NC Logical. Whether or not to include "Not.Correlated"
#' module. Defaults to \code{TRUE}.
#'
#' @return Module names
#' @examples
#' data("fc")
#' mod_names(fc)
#' @rdname mod_names
#' @export
setGeneric("mod_names", function(fc, include_NC = TRUE) {
  standardGeneric("mod_names")
})

#' @rdname mod_names
setMethod("mod_names", signature(fc = "fcoex"),
          function(fc, include_NC = TRUE) {
            mods <- NULL
            if (length(fc@module_list) > 0) {
              mods <- names(fc@module_list)
            } else {
              warning("No modules in this fcoex object.")
            }
            return(mods)
          })


#' Get the module genes in a fcoex object
#'
#' This function was copied and adapted from the CEMiTool package.
#'
#' @param fc Object of class \code{fcoex}
#' @param module A character string with the name of the module of which
#' genes are to be returned. Defaults to \code{NULL}, which returns the full
#' list of genes and modules.
#'
#' @return Object of class \code{data.frame} containing genes and their
#' respective module
#'
#' @rdname module_genes
#' @examples
#' data("fc")
#' module_genes(fc)
#' @export
setGeneric("module_genes", function(fc, module = NULL) {
  standardGeneric("module_genes")
})

#' @rdname module_genes
setMethod("module_genes", signature(fc = "fcoex"),
          function(fc, module = NULL) {
            # mod_names <- unique(fc@module[, "modules"])
            res <- NULL
            if (length(fc@module_list) > 0) {
              res <- fc@module_list
            } else {
              message("No modules in this fcoex object.")
              return(res)
            }
            mod_names <- names(fc@module_list)
            if (!is.null(module)) {
              if (module %in% mod_names) {
                res <- res[names(res) == module]
              } else {
                stop("Undefined module!")
              }
            }
            return(res)
          })

#' Print a fcoex object
#'
#' @param object Object of class fcoex
#'
#' @return A fcoex object.
#' @examples
#' data("fc")
#' fc
#' @export
setMethod("show", signature(object = "fcoex"),
          function(object) {
            cat("fcoex Object\n")
            cat("- Number of modules:", suppressWarnings(nmodules(object)),
                "\n")
            cat("- Module headers: \n")
            if (length(object@module_list) == 0) {
              cat("null\n")
            } else {
              cat(names(object@module_list), sep = ", ")
              cat("\n")
            }
            cat("- Expression file: ")
            if (nrow(object@expression) == 0) {
              cat("null\n")
            } else {
              cat(
                "data.frame with",
                nrow(object@expression),
                "genes and",
                ncol(object@expression),
                "cells\n"
              )
            }
            if (is.character(object@selected_genes)) {
              if (length(object@selected_genes) != nrow(object@expression)) {
                cat("- Selected data:",
                    length(object@selected_genes),
                    "genes selected\n")
              }
            }
          })

#' Retrieve over representation analysis (ORA) results
#'
#' @param fc Object of class \code{fcoex}
#'
#' @details This function returns the results of the \code{mod_ora} function
#' on the \code{fcoex} object. The ID column corresponds to pathways in
#' the gmt file for which  genes in the modules were enriched. The Count
#' column shows the number of genes in the module that are enriched for each
#' pathway. The GeneRatio column shows the proportion of  genes in the module
#' enriched for a given pathway out of all the genes in the module enriched
#' for any given pathway. The BgRatio column shows the proportion of genes
#' in a given pathway out of all the genes in the gmt file. For more details,
#' please refer to the \code{clusterProfiler} package documentation.
#'
#' This function was ipsis litteris adapted from the CEMiTool package.
#'
#' @return Object of class \code{data.frame} with ORA data
#' @examples
#' data("fc")
#' ora_data(fc)
#' @references Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler:
#' an R package for comparing biological themes among gene clusters. OMICS:
#' A Journal of Integrative Biology. 2012, 16(5):284-287.
#' @rdname ora_data
#' @export
setGeneric("ora_data", function(fc) {
  standardGeneric("ora_data")
})

#' @rdname ora_data
setMethod("ora_data", signature("fcoex"),
          function(fc) {
            return(fc@ora)
          })



#' Recluster cells based on fcoex module composition
#'
#' @param fc Object of class \code{fcoex}
#' @param hclust_method method for the hclust function. Defaults to "ward.D2".
#' @param dist_method  method for the dist function. Defaults to "manhattan".
#' @param k desired number of clustes. Defaults to 2.
#' @param verbose Adds verbosity, defaults to TRUE.
#' @return Object of class \code{data.frame} with new clusters
#' @examples
#' data("fc")
#' fc <- recluster(fc)
#' @export
#' @rdname recluster
setGeneric("recluster", function(fc,
                                 hclust_method = "ward.D2",
                                 dist_method = "manhattan",
                                 k = 2,
                                 verbose = TRUE) {
  standardGeneric("recluster")
})
#' @rdname recluster
setMethod("recluster", signature("fcoex"),
          function(fc,
                   hclust_method = "ward.D2",
                   dist_method = "manhattan",
                   k = 2,
                   verbose = TRUE) {
            mod_idents <- list()
            
            if (verbose) {
              message("Detecting clusters for the following modules: ")
            }
            
            for (i in names(fc@module_list)) {
              if (verbose) {
                message(print(i))
              }
              expression_table <-
                fc@expression[fc@module_list[[i]],]
              d <-
                dist(t(as.matrix(expression_table)), method = dist_method)
              hc <- hclust(d, method = hclust_method)
              idents <- as.factor(cutree(hc, k))
              
              if (k == 2) {
                mean_1 <- mean(as.numeric(fc@expression[i,][idents == 1]))
                mean_2 <-
                  mean(as.numeric(fc@expression[i,][idents == 2]))
                if (mean_1 > mean_2) {
                  # The first cluster will be the header positive cluster
                  first <- "HP"
                  second <- "HN"
                } else {
                  # The first cluster will be the header negative cluster
                  first <- "HN"
                  second <- "HP"
                }
                idents <- ifelse(idents == 1, first, second)
              }
              mod_idents[[i]] <- as.factor(idents)
            }
            fc@mod_idents <- mod_idents
            return(fc)
          })

#' Retrieves module identities from the recluster function
#'
#' @param fc Object of class \code{fcoex}
#'
#' @return Named object of class \code{list} with clusterings derived
#' from the recluster function.
#'
#' @examples
#' data("fc")
#' idents(fc)
#' @references Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler:
#' an R package for comparing biological themes among gene clusters. OMICS:
#' A Journal of Integrative Biology. 2012, 16(5):284-287.
#' @rdname idents
#' @export
setGeneric("idents", function(fc) {
  standardGeneric("idents")
})

#' @rdname idents
setMethod("idents", signature("fcoex"),
          function(fc) {
            return(fc@mod_idents)
          })

# Auxiliary functions ------

change_dots_for_dashes <- function(vector_of_genes) {
  gsub("\\.", "-", vector_of_genes)
}


get_gene_by_gene_correlation_matrix_in_series <-
  function(genes_from_su_ranking,
           expression_table_only_with_genes_with_high_su) {
    gene_by_gene_su_correlation <-
      data.frame(genes = genes_from_su_ranking)
    
    pb_findclusters <-
      progress_bar$new(total = length(genes_from_su_ranking),
                       format = "[:bar] :percent eta: :eta")
    for (gene_i in genes_from_su_ranking) {
      pb_findclusters$tick()
      
      discrete_vector_of_gene_i <-
        as.factor(expression_table_only_with_genes_with_high_su[gene_i,])
      
      gene_i_correlates <-
        FCBF::get_su_for_feature_table_and_vector(
          feature_table = expression_table_only_with_genes_with_high_su,
          target_vector = as.factor(discrete_vector_of_gene_i)
        )
      
      gene_i_correlates$gene <-
        change_dots_for_dashes(rownames(gene_i_correlates))
      
      # Reordering rows
      gene_i_correlates <-
        gene_i_correlates[match(gene_by_gene_su_correlation$genes,
                                gene_i_correlates$gene),]
      
      colnames(gene_i_correlates)[1] <- gene_i
      
      gene_by_gene_su_correlation[, gene_i] <-
        gene_i_correlates[, 1]
    }
    
    gene_by_gene_su_correlation <- gene_by_gene_su_correlation[,-1]
    
    
    return(gene_by_gene_su_correlation)
  }

get_gene_by_gene_correlation_matrix_in_parallel <-
  function(genes_from_su_ranking,
           expression_table_only_with_genes_with_high_su,
           discretized_exprs,
           su_to_class_higher_than_minimum_su) {
    # get and adjacency matrix for gene to gene correlation
    gene_by_gene_su_correlation <-
      data.frame(genes = genes_from_su_ranking)
    
    cl <- detectCores() - 2
    bla <- mclapply(genes_from_su_ranking, function(i) {
      .get_correlates(
        i,
        gene_by_gene_su_correlation,
        discretized_exprs,
        expression_table_only_with_genes_with_high_su
      )
    }, mc.cores = cl)
    gene_by_gene_su_correlation <- as.data.frame(bla)
    rownames(gene_by_gene_su_correlation) <-
      su_to_class_higher_than_minimum_su$gene
    colnames(gene_by_gene_su_correlation) <-
      su_to_class_higher_than_minimum_su$gene
    
    return(gene_by_gene_su_correlation)
  }

trim_correlation_matrix <- function(genes_from_su_ranking,
                                    gene_by_gene_su_correlation,
                                    su_to_class,
                                    su_to_class_higher_than_minimum_su) {
  filtered_gene_by_gene_su_correlation <-
    data.frame(genes = genes_from_su_ranking)
  
  for (i in colnames(gene_by_gene_su_correlation)) {
    tf_vector <-
      gene_by_gene_su_correlation[, i] > su_to_class$SU[seq_along(su_to_class_higher_than_minimum_su$gene)]
    filtered_gene_by_gene_su_correlation[, i] <-
      gene_by_gene_su_correlation[, i] * tf_vector
  }
  
  return(filtered_gene_by_gene_su_correlation)
}


get_list_of_modules <-
  function(genes_from_fcbf_filter,
           filtered_gene_by_gene_su_correlation) {
    list_of_fcbf_modules <- list()
    
    for (seed in genes_from_fcbf_filter) {
      
      correlated_genes_score <-
        filtered_gene_by_gene_su_correlation[, seed]
      
      module_members <-
        as.character(filtered_gene_by_gene_su_correlation$genes[correlated_genes_score > 0])
      
      module_members <- module_members[!is.na(module_members)]
      if (length(module_members) > 1) {
        list_of_fcbf_modules[[seed]] <- module_members
      }
    }
    
    return(list_of_fcbf_modules)
  }


run_fcbf_for_module_headers <- function(discretized_exprs,
                                        target,
                                        n_genes_selected_in_first_step,
                                        minimum_su,
                                        verbose = verbose) {
  output_of_fcbf_filter <-
    FCBF::fcbf(
      discretized_exprs,
      target,
      n_genes_selected_in_first_step,
      minimum_su = minimum_su,
      verbose = verbose
    )
  
  output_of_fcbf_filter$gene <- rownames(output_of_fcbf_filter)
  
  output_of_fcbf_filter$gene <-
    change_dots_for_dashes(output_of_fcbf_filter$gene)
  
  genes_from_fcbf_filter <- output_of_fcbf_filter$gene
  
  return(genes_from_fcbf_filter)
}
