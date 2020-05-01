# Many parts of this code were directly adapted from the
# CEMiTool package. This include chunks of copied and pasted code
# located inside the source code for CEMiTool functions. 
# Functions that contained adapted code were explicit denoted.



#' @importFrom grDevices rainbow dev.off pdf
#' @importFrom utils write.table head
#' @importFrom stats cutree dist hclust
#' @importFrom methods new 'slot<-' show
#' @importFrom pathwayPCA read_gmt
#' @import SingleCellExperiment
#' @import dplyr


setOldClass('gg')
setOldClass('ggplot')
setOldClass('gtable')

#' An S4 class to represent the fcoex analysis.
#'
#' @slot expression Normalized gene expression table from 
#' single-cells \code{data.frame}.
#' @slot discretized_expression Discretized gene expression table from 
#' single-cells \code{data.frame}.
#' @slot target Original target classes for the cells (\code{factor}).
#' @slot selected_genes Character \code{vector} containing the names of 
#' genes selected for analysis
#' @slot module_list \code{list} containing genes in each module.
#' @slot adjacency \code{data.frame} containing the adjacency table for 
#' the selected genes before trimming.
#' @slot adjacency_trimmed \code{data.frame} containing the adjacency table for 
#' the selected genes after trimming. 
#' @slot coex_network_plot list of ggplot graphs with module gene interactions.
#' @slot new_clusters \code{list} containing gene interactions present in 
#' modules.
#' @slot mod_colors character \code{vector} containing colors associated with 
#' each network module.
#' @slot ora Over-representation analysis results \code{data.frame}.
#' @slot barplot_ora list of ggplot graphs with over-representation analysis 
#' results per module.
#' @slot mod_idents Identities of cells based on each co-expression module. 
#' Determined by the "recluster" method
#' @slot parameters \code{list} containing analysis parameters.
setClass(
  'fcoex',
  slots = list(
    expression = 'data.frame',
    discretized_expression = 'data.frame',
    target = 'factor',
    selected_genes = 'vector',
    module_list = 'list',
    adjacency = 'list',
    adjacency_trimmed = 'list',
    coex_network_plot = 'list',
    new_clusters = 'list',
    mod_colors = 'character',
    parameters = 'list',
    ora = 'data.frame',
    barplot_ora = 'list',
    mod_idents = 'list'
  )
)

setMethod("initialize", signature = "fcoex",
          function(.Object, expression, target) {
            .Object@expression <- expression
            .Object@target <- target
            return(.Object)
          })

#' Create a fcoex object
#'
#' @param expression Normalized gene expression table from single-cells
#'  \code{data.frame}.
#' @param target Original target classes for the cells (\code{factor}).
#' @return Object of class \code{fcoex}
#' @examples
#' # Create new fcoex object
#' library(SingleCellExperiment)
#' data("mini_pbmc3k")
#' targets <- colData(mini_pbmc3k)$clusters
#' exprs <- as.data.frame(assay(mini_pbmc3k, "logcounts"))
#' fc <- new_fcoex(exprs, targets)

#' @export
new_fcoex <- function(expression = data.frame(), target = vector()) {
  fc <- new("fcoex", expression = expression, target = target)
  msg <- "Created new fcoex object."
  message(msg)
  return(fc)
}

#' Set the discretized expression attribute
#' Uses the discretize_exprs function of the FCBF package
#'
#' @param fc Object of class \code{fcoex}
#' @param method Method applied to all genes for discretization. 
#' Methods available: "varying_width"
#'  (Binarization modulated by the number_of_bins param),
#' "mean" (Split in ON/OFF by each gene mean expression),
#' "median" (Split in ON/OFF by each gene median expression),
#' "mean_sd"(Split in low/medium/high by each assigning "medium" to 
#' the interval between mean +- standard_deviation.
#' Modulated by the alpha param, which enlarges (>1) or shrinks (<1) the 
#' "medium" interval. ),
#' ),
#' "kmeans"(Split in different groups by the kmeans algorithm. As many 
#' groups as specified by the centers param) and
#' "min_max_\%" (Similat to the "varying width", a binarization threshold 
#' in a % of the min-max range is set. (minmax\% param)),
#' "GMM" (A Gaussian Mixture Model as implemented by the package mclust, 
#' trying to fit 2:5 Gaussians). Default is "varying_width"
#' 
#' @param number_of_bins Number of equal-width bins for discretization.
#' Note: it is a binary discretization, with the
#' first bin becoming one class ('low') and the other bins, another class 
#' ('high').#' Defaults to 4.
#' @param alpha Modulator for the "mean_sd" method.Enlarges (>1) or 
#' shrinks (<1) the "medium" interval. Defaults to 1.
#' @param centers Modulator for the "kmeans" method. Defaults to 3.
#' @param min_max_cutoff <- Modulator for the "min_max_\%" method. 
#' Defaults to 0.25.
#' @param show_pb Enables a progress bar for the discretization. 
#' Defaults to TRUE.
#' @return A data frame with the discretized features in the same 
#' order as previously
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
setGeneric("discretize", function(fc, number_of_bins = 4,
                                  method = "varying_width",
                                  alpha = 1,
                                  centers = 3,
                                  min_max_cutoff = 0.25,
                                  show_pb = TRUE) {
  standardGeneric("discretize")
})

#' @rdname discretize
setMethod("discretize", signature("fcoex"),
          function(fc,
                   number_of_bins = 4,
                   method = "varying_width",
                   alpha = 1,
                   centers = 3,
                   min_max_cutoff = 0.25,
                   show_pb = TRUE) {
            expression_table <- fc@expression
            discretized_expression <-
              FCBF::discretize_exprs(expression_table,
                                     number_of_bins,
                                     method,
                                     alpha,
                                     centers,
                                     min_max_cutoff,
                                     show_pb)
            colnames(discretized_expression) <-
              colnames(expression_table)
            fc@discretized_expression <- discretized_expression
            
            return(fc)
          })


#' .get_correlates 
#' 
#' auxiliary function for find_cbf_modules
#' @param i A gene to be correlated
#' @param su_i_j_matrix the dataframe with the correlations to be updated
#' @param discretized_exprs the dataframe with discretized expression 
#' to extract a gene
#' @param exprs_small the dataframe to after the filtering step
#' @return the updated column of the su_i_j_matrix
.get_correlates <- function(i,
           su_i_j_matrix,
           discretized_exprs,
           exprs_small) {
    gene_i <- as.factor(discretized_exprs[i,])
    gene_i_correlates <-
      FCBF::get_su_for_feature_table_and_vector(feature_table = exprs_small, target_vector = as.factor(exprs_small[i,]))
    gene_i_correlates$gene <-
      gsub('\\.', '-', rownames(gene_i_correlates))
    gene_i_correlates <-
      gene_i_correlates[match(su_i_j_matrix$genes, gene_i_correlates$gene), ]
    colnames(gene_i_correlates)[1] <- i
    su_i_j_matrix[, i] <- gene_i_correlates[, 1]
    su_i_j_matrix[, i]
  }

#' find_cbf_modules
#'
#' find_cbf_modules uses Symmetrical Uncertainty as a correlation measure
#'  and the FCBF algorithm to
#'
#' 1 - Filter the gene list by correlations to a class (Step 1)
#'
#' and
#'
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
  
  discretized_exprs <- fc@discretized_expression
  target <- fc@target
  
  # get the SU scores for each gene
  message('Getting SU scores')
  su_to_class <- FCBF::get_su_for_feature_table_and_vector(discretized_exprs, target)
  su_to_class$gene <- gsub('\\.', '-', rownames(su_to_class))
  colnames(su_to_class)[1] <- 'SU'
  
  
  # Run FCBF with the parameters of the function. 
  message('Running FCBF to find module headers')
  fcbf_filtered <-
    FCBF::fcbf(discretized_exprs,
               target,
               n_genes_selected_in_first_step,
               minimum_su = FCBF_threshold,
               verbose = verbose)
  
  fcbf_filtered$gene <- rownames(fcbf_filtered)
  # R does not like points. Subs for -.
  FCBF_genes <- gsub('\\.', '-', fcbf_filtered$gene)
  if (length(n_genes_selected_in_first_step)) {
    FCBF_threshold <- su_to_class$SU[n_genes_selected_in_first_step]
  }
  
  # get only those with an SU score above a threshold
  SU_threshold <- FCBF_threshold
  su_to_class_small <-
    su_to_class[su_to_class[1] > SU_threshold, ]
  
  SU_genes <- gsub('\\.', '-', su_to_class_small[, 2])
  
  # Save the selected SU_genes  in the fc object
  fc@selected_genes <- SU_genes
  
  exprs_small <- discretized_exprs[SU_genes , ]
  

  # get and adjacency matrix for gene to gene correlation
  su_i_j_matrix <- data.frame(genes =  SU_genes)
  message('Calculating adjacency matrix')
  
  if (!is_parallel) {
    pb_findclusters <- progress_bar$new(total = length(SU_genes),
                                        format =   "[:bar] :percent eta: 
                                        :eta")
    # this can surely be improved for speed.
    for (i in SU_genes) {
      pb_findclusters$tick()
      gene_i <- as.factor(discretized_exprs[i,])
      gene_i_correlates <-
        FCBF::get_su_for_feature_table_and_vector(feature_table = exprs_small, target_vector = as.factor(exprs_small[i,]))
      gene_i_correlates$gene <-
        gsub('\\.', '-', rownames(gene_i_correlates))
      gene_i_correlates <-
        gene_i_correlates[match(su_i_j_matrix$genes, gene_i_correlates$gene), ]
      colnames(gene_i_correlates)[1] <- i
      su_i_j_matrix[, i] <- gene_i_correlates[, 1]
      
    }
    su_i_j_matrix <- su_i_j_matrix[, -1]

  }  else {
    cl <- detectCores() - 2
    bla <- mclapply(SU_genes, function(i) {
      .get_correlates(i, su_i_j_matrix, discretized_exprs, exprs_small)
    }, mc.cores = cl)
    su_i_j_matrix <- as.data.frame(bla)
    rownames(su_i_j_matrix) <- su_to_class_small$gene
    colnames(su_i_j_matrix) <- su_to_class_small$gene
  }
  
  filtered_su_i_j_matrix <- data.frame(genes =  SU_genes)
  
  message('Trimming and getting modules from adjacency matrix')
  for (i in colnames(su_i_j_matrix)) {
    tf_vector <-
      su_i_j_matrix[, i] > su_to_class$SU[seq_along(su_to_class_small$gene)]
    filtered_su_i_j_matrix[, i] <- su_i_j_matrix[, i] * tf_vector
  }
  
  
  list_of_fcbf_modules <- list()
  for (seed in FCBF_genes) {
    module_members <-
      as.character(filtered_su_i_j_matrix$genes[filtered_su_i_j_matrix[, seed] >
                                                  0])
    module_members <- module_members[!is.na(module_members)]
    if (length(module_members) > 1) {
      list_of_fcbf_modules[[seed]] <- module_members
    }
  }
  
  fc@adjacency <- su_i_j_matrix
  fc@adjacency_trimmed <- filtered_su_i_j_matrix
  fc@module_list <- list_of_fcbf_modules
  return(fc)
})



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
#'
#' @export
setGeneric('nmodules', function(fc) {
  standardGeneric('nmodules')
})

#' @rdname nmodules
setMethod('nmodules', signature('fcoex'),
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
setGeneric('mod_gene_num', function(fc, module = NULL) {
  standardGeneric('mod_gene_num')
})
#' @rdname mod_gene_num
setMethod('mod_gene_num', signature(fc = 'fcoex'),
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
setGeneric('mod_names', function(fc, include_NC = TRUE) {
  standardGeneric('mod_names')
})

#' @rdname mod_names
setMethod('mod_names', signature(fc = 'fcoex'),
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
setGeneric('module_genes', function(fc, module = NULL) {
  standardGeneric('module_genes')
})

#' @rdname module_genes
setMethod('module_genes', signature(fc = 'fcoex'),
          function(fc, module = NULL) {
            #mod_names <- unique(fc@module[, "modules"])
            res <- NULL
            if (length(fc@module_list) > 0) {
              res <- fc@module_list
            } else{
              message("No modules in this fcoex object.")
              return(res)
            }
            mod_names <- names(fc@module_list)
            if (!is.null(module)) {
              if (module %in% mod_names) {
                res <- res[names(res) == module]
              } else{
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
setMethod('show', signature(object = 'fcoex'),
          function(object) {
            cat("fcoex Object\n")
            cat("- Number of modules:", suppressWarnings(nmodules(object)), 
                "\n")
            cat("- Module headers: \n")
            if (length(object@module_list) == 0) {
              cat("null\n")
            } else {
              cat(names(object@module_list), sep = ", ")
              cat('\n')
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

#' # Run module overrepresentation analysis
#'
#' This function was modified from the CEMiTool package.
#' Chunks of code were retained "as is"
#' 
#' @param fc A fcoex object.
#' @param gmt A gmt file with gene sets for ora analysis
#' @param verbose Controls verbosity. Defaults to FALSE. 
#' @examples 
#' data("fc")
#' gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
#' gmt_in <- pathwayPCA::read_gmt(gmt_fname)
#' fc <- mod_ora(fc, gmt_in)
#' @return  A fcoex object containing over-representation analysis data
#' @rdname mod_ora
#' @export
setGeneric('mod_ora', function(fc, gmt, verbose = FALSE) {
  standardGeneric('mod_ora')
})

#' @rdname mod_ora
setMethod('mod_ora', signature('fcoex'),
          function(fc, gmt, verbose = FALSE) {
            #fc <- get_args(fc, vars=mget(ls()))
            if (!is(gmt, "pathwayCollection")) {
              stop("The gmt object must be loaded via pathwayPCA::read_gmt and be a pathwayCollectionobject")
            }
          
            pathwayPCA_gmt <- gmt
            gmt_df <- pathwayPCA_gmt$pathways

            gmt_df <- as.data.frame(unlist(gmt_df), use.names = TRUE)
            pathway_sizes <- unlist(lapply(pathwayPCA_gmt$pathways, length))
            colnames(gmt_df) <- "gene"
            gmt_df$term <- rep(pathwayPCA_gmt$TERMS, pathway_sizes)

            
            gmt_df<-gmt_df[,c("term","gene")]
            
            
            if (verbose) {
              message('Running ORA')
              message("Using all genes in GMT file as universe.")
            }
            allgenes <- unique(gmt_df[, "gene"])
            if (is.null(fc@module_list)) {
              warning("No modules in fcoex object! Did you run find_cbf_modules()?")
              return(fc)
            }
            mods <- fc@module_list
            res_list <-
              lapply(names(mods), ora, gmt_df, allgenes, mods)
            if (all(lapply(res_list, nrow) == 0)) {
              warning(
                "Enrichment is NULL. Either your gmt file is inadequate or your modules really aren't enriched for any of the pathways in the gmt file."
              )
              return(fc)
            }
            names(res_list) <- names(mods)
            
            res <- lapply(names(res_list), function(x) {
              if (nrow(res_list[[x]]) > 0) {
                as.data.frame(cbind(x, res_list[[x]]))
              }
            })
            res <- do.call(rbind, res)
            names(res)[names(res) == "x"] <- "Module"
            
            rownames(res) <- NULL
            fc@ora <- res
            return(fc)
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
setGeneric("recluster", function(fc, hclust_method = "ward.D2",
                                 dist_method = 'manhattan',
                                 k = 2,
                                 verbose = TRUE) {
  standardGeneric("recluster")
})
#' @rdname recluster
setMethod("recluster", signature("fcoex"),
          function(fc,
                   hclust_method = "ward.D2",
                   dist_method = 'manhattan',
                   k = 2,
                   verbose = TRUE) {
            mod_idents <- list()
  
            if (verbose){
              message("Detecting clusters for the following modules: ")   
            }

            for (i in names(fc@module_list)) {
              if (verbose){
                message(print(i))   
              }
              expression_table <-
                fc@expression[fc@module_list[[i]], ]
              d <-
                dist(t(as.matrix(expression_table)), method = dist_method)
              hc <- hclust(d, method = hclust_method)
              idents <- as.factor(cutree(hc, k))

              if (k == 2){
               mean_1 <-  mean(as.numeric(fc@expression[i,][idents==1]))
               mean_2 <-  mean(as.numeric(fc@expression[i,][idents==2]))
               if (mean_1 > mean_2){
                 # The first cluster will be the header positive cluster
                 first = "HP"
                 second = "HN"
               } else {
                 # The first cluster will be the header negative cluster
                 first = "HN"
                 second = "HP"
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
