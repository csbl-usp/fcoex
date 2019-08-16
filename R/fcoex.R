#' @importFrom grDevices rainbow
#' @importFrom utils write.table
#' @importFrom methods new 'slot<-' show
#' @import dplyr

setOldClass('gg')
setOldClass('ggplot')
setOldClass('gtable')

#' An S4 class to represent the fcoex analysis.
#'
#' @slot expression Normalized gene expression table from single-cells \code{data.frame}.
#' @slot discretized_expression Discretized gene expression table from single-cells \code{data.frame}.
#' @slot target Original target classes for the cells (\code{factor}).
#' @slot selected_genes Character \code{vector} containing the names of genes  selected for analysis
#' @slot module_list \code{list} containing genes in each module.
#' @slot adjacency \code{data.frame} containing the adjacency table for the selected genes.
#' @slot interaction_plot list of ggplot graphs with module gene interactions.
#' @slot new_clusters \code{list} containing gene interactions present in modules.
#' @slot mod_colors character \code{vector} containing colors associated with each network module.
#' @slot ora Over-representation analysis results \code{data.frame}.
#' @slot barplot_ora list of ggplot graphs with over-representation analysis results per module.
#' @slot mod_idents Identities of cells based on each co-expression module. Determined by the "recluster" method
# #' @slot parameters \code{list} containing analysis parameters.
#' @examples
#' # Get example expression data
#' data(expr0)
#' # Initialize fcoex object with expression
#' fc <- new("fcoex", expression=expr0)
setClass('fcoex', slots=list(expression='data.frame',
                             discretized_expression ='data.frame',
                                target='factor',
                                selected_genes='vector',
                                module_list='list',
                                adjacency='list',
                                interaction_plot='list',
                                new_clusters='list',
                                mod_colors='character',
                                parameters='list',
                                ora='data.frame',
                                barplot_ora='list',
                                mod_idents='list'))

setMethod("initialize", signature="fcoex",
    function(.Object, expression,target){
        .Object@expression <- expression
        .Object@target <- target
        return(.Object)
    })

#' Create a fcoex object
#'
#' @param expression Normalized gene expression table from single-cells \code{data.frame}.
#' @param target Original target classes for the cells (\code{factor}).
#' @return Object of class \code{fcoex}
#' @examples
#' # Create new fcoex object
#' fc <- new_fcoex()
#' @export
new_fcoex <- function(expr=data.frame(), target=vector()){
    fc <- new("fcoex", expression=expr, target=target)
    msg <- "Created new fcoex object."
    message(msg)
    return(fc)
}

#' Set the discretized expression attribute
#' Uses the discretize_exprs function of the FCBF package
#'
#' @param fc Object of class \code{fcoex}
#' @param method Method applied to all genes for discretization. Methods available: "varying_width"
#'  (Binarization modulated by the number_of_bins param),
#' "mean" (Split in ON/OFF by each gene mean expression),
#' "median" (Split in ON/OFF by each gene median expression),
#' "mean_sd"(Split in low/medium/high by each assigning "medium" to the interval between mean +- standard_deviation.
#' Modulated by the alpha param, which enlarges (>1) or shrinks (<1) the "medium" interval. ),
#' ),
#' "kmeans"(Split in different groups by the kmeans algorithm. As many groups as specified by the centers param) and
#' "min_max_\%" (Similat to the "varying width", a binarization threshold in a % of the min-max range is set. (minmax\% param)),
#' "GMM" (A Gaussian Mixture Model as implemented by the package mclust, trying to fit 2:5 Gaussians). Default is "varying_width"
#' @param number_of_bins Number of equal-width bins for discretization.
#' Note: it is a binary discretization, with the
#' first bin becoming one class ('low') and the other bins, another class ('high').
#' Defaults to 4.
#' @param alpha Modulator for the "mean_sd" method.Enlarges (>1) or shrinks (<1) the "medium" interval. Defaults to 1.
#' @param centers Modulator for the "kmeans" method. Defaults to 3.
#' @param min_max_cutoff <- Modulator for the "min_max_\%" method. Defaults to 0.25.
#' @param show_pb Enables a progress bar for the discretization. Defaults to TRUE.
#' @return A data frame with the discretized features in the same order as previously
#' @import FCBF
#' @export
#' @rdname discretize
setGeneric("discretize", function(fc, ...) {
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
                   show_pb = TRUE){
             expression_table <- fc@expression
             discretized_expression <-FCBF::discretize_exprs(expression_table,
                                                            number_of_bins,
                                                            method,
                                                            alpha,
                                                            centers,
                                                            min_max_cutoff,
                                                            show_pb)
            colnames(discretized_expression) <- colnames(expression_table)
            fc@discretized_expression <- discretized_expression
            
            return(fc)
          })

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
#' @param is_parallel Uses package parallel to paralleliza calculations. Defaults to FALSE.
#' @param verbose Adds verbosity. Defaults to TRUE
#' @param n_genes Sets the number of genes to be selected in the first part of the algorithm.
#' If left unchanged, it defaults to NULL and the thresh parameter is used.
#' Caution: it overrides the thresh parameter altogether.
#' @return Returns a list with the CBF modules found or a adjacency matrix of the graph
#' @import dplyr
#' @import parallel
#' @import progress
#' @import FCBF
#' @export
#' @rdname find_cbf_modules
setGeneric("find_cbf_modules", function(fc, ...) {
standardGeneric("find_cbf_modules")
})

#' @rdname find_cbf_modules
setMethod("find_cbf_modules", signature("fcoex"), function(fc, 
                                                           n_genes = NULL, 
                                                           FCBF_threshold = 0.1, 
                                                           verbose = TRUE,
                                                           is_parallel = FALSE){
  discretized_exprs <- fc@discretized_expression
  target <- fc@target
  
  # get the SU scores for each gene
  message('Getting SU scores')
  su_ic_vector <- FCBF::get_su(discretized_exprs, target)
  su_ic_vector$gene <- gsub('\\.', '-',rownames(su_ic_vector))
  
  colnames(su_ic_vector)[1] <- 'SU'
  message('Running FCBF to find module headers')
  fcbf_filtered <- FCBF::fcbf(discretized_exprs, target, n_genes, thresh = FCBF_threshold, verbose = verbose)
  fcbf_filtered$gene <- rownames(fcbf_filtered)
  # R does not like points. Subs for -.
  FCBF_genes <- gsub('\\.', '-', fcbf_filtered$gene)
  
  if (length(n_genes)){
    FCBF_threshold <- su_ic_vector$SU[n_genes]
  }
  # get only those with an SU score above a threshold
  SU_threshold <- FCBF_threshold
  
  su_ic_vector_small <- su_ic_vector[su_ic_vector[1] > SU_threshold,]
  
  
  SU_genes <- gsub('\\.', '-',su_ic_vector_small[,2])
  
  fc@selected_genes <-SU_genes
  
  exprs_small <- discretized_exprs[SU_genes ,]
  
  
  # get and adjacency matrix for gene to gene correlation
  su_i_j_matrix <- data.frame(genes =  SU_genes)
  message('Calculating adjacency matrix')
  
  if (!is_parallel){
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
  su_i_j_matrix <- su_i_j_matrix[,-1]
  }
  
  #      This was not faster than the for loop! ######
  #  
#      get_correlates <- function(i, su_i_j_matrix, discretized_exprs, exprs_small){
#    gene_i <- as.factor(discretized_exprs[i, ])
#     gene_i_correlates <- FCBF::get_su(x = exprs_small, y = as.factor(exprs_small[i, ]))
#     gene_i_correlates$gene <- gsub('\\.', '-',rownames(gene_i_correlates))
#   gene_i_correlates <- gene_i_correlates[match(su_i_j_matrix$genes,gene_i_correlates$gene),]
#     colnames(gene_i_correlates)[1] <- i
#      su_i_j_matrix[, i] <- gene_i_correlates[,1]
#      su_i_j_matrix[, i]
#    }
  #   
  #   bla <- pblapply(SU_genes, function(x){
  #     get_correlates(x, su_i_j_matrix, discretized_exprs, exprs_small)
  #   })
  #   
  #   su_i_j_matrix <- as.data.frame(ble)
  ################################################## 
  
  # Second Try
  else{
  get_correlates <- function(i, su_i_j_matrix, discretized_exprs, exprs_small){
    gene_i <- as.factor(discretized_exprs[i, ])
    gene_i_correlates <- FCBF::get_su(x = exprs_small, y = as.factor(exprs_small[i, ]))
    gene_i_correlates$gene <- gsub('\\.', '-',rownames(gene_i_correlates))
    gene_i_correlates <- gene_i_correlates[match(su_i_j_matrix$genes,gene_i_correlates$gene),]
    colnames(gene_i_correlates)[1] <- i
    su_i_j_matrix[, i] <- gene_i_correlates[,1]
    su_i_j_matrix[, i]
  }
  cl <- detectCores()-2
  bla <- mclapply(SU_genes, function(i){
    get_correlates(i, su_i_j_matrix, discretized_exprs, exprs_small)
 }, mc.cores = cl)
  su_i_j_matrix <- as.data.frame(bla)
  rownames(su_i_j_matrix) <- su_ic_vector_small$gene
  colnames(su_i_j_matrix) <- su_ic_vector_small$gene
  }
  
  filtered_su_i_j_matrix <- data.frame(genes =  SU_genes)
  
  message('Getting modules from adjacency matrix')
  for (i in colnames(su_i_j_matrix)){
      tf_vector <- su_i_j_matrix[,i] > su_ic_vector$SU[seq_len(length(su_ic_vector_small$gene))]
  
      filtered_su_i_j_matrix[,i] <- su_i_j_matrix[,i] * tf_vector
  }
  
  
  list_of_fcbf_modules <- list()
  for (seed in FCBF_genes){
    module_members <- as.character(filtered_su_i_j_matrix$genes[filtered_su_i_j_matrix[,seed]>0])
    module_members <- module_members[!is.na(module_members)]
    if(length(module_members) > 1) {
      list_of_fcbf_modules[[seed]] <- module_members
    }
  }
  
  fc@adjacency <- filtered_su_i_j_matrix
  fc@module_list <- list_of_fcbf_modules
  return(fc)
})



#' Set module colors mod_colors attribute
#' @param fc Object of class \code{fcoex}
#'
#' @return A vector with color names.
#' @rdname mod_colors
#' @export
setGeneric("mod_colors", function(fc) {
    standardGeneric("mod_colors")
})

#' @rdname mod_colors
setMethod("mod_colors", signature("fcoex"),
    function(fc){
       mod_names <- names(fc@module_list)
       nmod <- length(mod_names)
       cols <- fc@mod_colors
       if(nmod != 0) {
           if(length(fc@mod_colors) == 0){
               if(nmod <= 16) {
                   cols <- rainbow(16, s = 1, v = 0.7)[1:nmod]
               } else {
                   cols <- rep(rainbow(16, s = 1, v = 0.7), ceiling(nmod/16))[1:nmod]
               }
               names(cols) <- mod_names
           } else {
               if(is.null(names(fc@mod_colors))){
                   warning("mod_colors should be a character vector with names corresponding to the modules")
               } else if(!all(sort(names(fc@mod_colors)) == sort(mod_names))){
                   warning("mod_colors names do not match with modules!")
               }
           }
       }
       fc@mod_colors <- cols
       return(fc)
    })

#' Get the number of modules in a fcoex object
#'
#' @param fc Object of class \code{fcoex}
#'
#' @return number of modules
#'
#' @rdname nmodules
#' @examples
#' # Get example fcoex object
#' data(fc)
#' # Get the number of modules
#' nmodules(fc)
#'
#' @export
setGeneric('nmodules', function(fc) {
    standardGeneric('nmodules')
})

#' @rdname nmodules
setMethod('nmodules', signature('fcoex'),
    function(fc){
        n <- 0
        if ((length(fc@module_list)) > 0){
            n <- length(fc@module_list)
        } else {
            warning("Run find_cbf_modules function to get modules!")
        }
        return(n)
    })




#' Get the number of genes in modules in a fcoex object
#'
#' @param fc Object of class \code{fcoex}
#' @param module Default is NULL. If a character string designating a module is
#' given, the number of genes in that module is returned instead.
#' @return The number of genes in module(s)
#'
#' @rdname mod_gene_num
#' @export
setGeneric('mod_gene_num', function(fc, module=NULL) {
    standardGeneric('mod_gene_num')
})
#' @rdname mod_gene_num
setMethod('mod_gene_num', signature(fc='fcoex'),
         function(fc, module=NULL){
              if(!is.null(module)){
                  if(!(all(module %in% mod_names(fc)))){
                      stop("Module '", module, "' not in fcoex object!")
                  }
              }
              if(!nrow(module_genes(fc))>0){
                  stop("No modules in fcoex object!")
              }
              if(!is.null(module)){
                  mod_genes <- fc@module_list[[module]]
              }
              return(mod_genes)
          })



#' Get module names in a fcoex object
#'
#' @param fc Object of class \code{fcoex}
#' @param include_NC Logical. Whether or not to include "Not.Correlated"
#' module. Defaults to \code{TRUE}.
#'
#' @return Module names
#'
#' @rdname mod_names
#' @examples
#' # Get example fcoex object
#' data(fc)
#' # Get module names
#' mod_names(fc)
#'
#' @export
setGeneric('mod_names', function(fc, include_NC=TRUE) {
    standardGeneric('mod_names')
})

#' @rdname mod_names
setMethod('mod_names', signature(fc='fcoex'),
          function(fc, include_NC=TRUE) {
              mods <- NULL
              if(length(fc@module_list) > 0){
                  mods <- names(fc@module_list)
              } else {
                  warning("No modules in this fcoex object.")
              }
              return(mods)
          }
)


#' Get the module genes in a fcoex object
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
#' # Get example fcoex object
#' data(fc)
#' # Get the module genes
#' module_genes(fc)
#' # Get genes for module M1
#' module_genes(fc, module="M1")
#' @export
setGeneric('module_genes', function(fc, module=NULL) {
    standardGeneric('module_genes')
})

#' @rdname module_genes
setMethod('module_genes', signature(fc='fcoex'),
    function(fc, module=NULL){
        #mod_names <- unique(fc@module[, "modules"])
        res <- NULL
        if(length(fc@module_list) > 0){
            res <- fc@module_list
        }else{
            message("No modules in this fcoex object.")
            return(res)
        }
        mod_names <- names(fc@module_list)
        if(!is.null(module)){
            if(module %in% mod_names){
                res <- res[names(res)==module]
            }else{
                stop("Undefined module!")
            }
        }
        return(res)
    }
)

#' Print a fcoex object
#'
#' @param object Object of class fcoex
#'
#' @return A fcoex object.
#'
#' @export
setMethod('show', signature(object='fcoex'),
    function(object) {
        cat("fcoex Object\n")
        cat("- Number of modules:", suppressWarnings(nmodules(object)), "\n")
        cat("- Module headers: \n")
        if(length(object@module_list) == 0){
            cat("null\n")
        } else {
            cat(names(fc@module_list), sep = ", ")
            cat('\n')
        }
        cat("- Expression file: ")
        if(nrow(object@expression) == 0){
            cat("null\n")
        } else {
            cat("data.frame with", nrow(object@expression), "genes and", ncol(object@expression), "cells\n")
        }
        if(is.character(object@selected_genes)){
            if(length(object@selected_genes) != nrow(object@expression)){
                cat("- Selected data:", length(object@selected_genes), "genes selected\n")
            }
        }
    }
)

#' Transform module genes list to a gmt file
#'
#' @keywords internal
#'
#' @param fc
#'
#' @return A .gmt file containing module genes in each row
#'
module_to_gmt <- function(fc, directory="./Tables"){
    if(length(fc@module_list) == 0){
        stop("No modules in fcoex object! Did you run find_modules()?")
    }else{
        gene_modules <- fc@module
        n_genes <- as.numeric(table(gene_modules$modules[gene_modules$modules != "Not.Correlated"]))
        n_genes <- n_genes[1:(length(n_genes))]
        module_names <- as.character(unique(gene_modules[, "modules"]))
        module_names <- module_names[module_names != "Not.Correlated"]
        module_names <- module_names[order(nchar(module_names), module_names)]

        gmt_df  <- as.data.frame(matrix("", ncol = max(n_genes), nrow = length(n_genes)), stringsAsFactors = FALSE)

        rownames(gmt_df) <- module_names

        for (i in 1:length(module_names)){
            mod <- module_names[i]
            selected <- gene_modules[gene_modules$modules == mod, "genes"]
                gmt_df[mod, 1:(length(selected))] <- selected
        }

        gmt_df <- as.data.frame(cbind(module_names, gmt_df))
        write.table(gmt_df, file.path(directory, "modules_genes.gmt"), sep="\t", col.names = FALSE, quote=FALSE)
    }
}


#' Save the fcoex object in files
#'
#' @param fc Object of class \code{fcoex}
#' @param directory a directory
#' @param force if the directory exists the execution will not stop
#' @param ... Optional parameters
#' @return A directory containing fcoex results in files.
#' @examples
#' # Get example fcoex object
#' data(fc)
#' # Save fcoex results in files
#' write_files(fc, directory=".", force=TRUE)
#'
#' @rdname write_files
#' @export
setGeneric('write_files', function(fc, ...) {
    standardGeneric('write_files')
})

#' @rdname write_files
setMethod('write_files', signature(fc='fcoex'),
    function(fc, directory="./Tables", force=FALSE) {
        if(dir.exists(directory)){
            if(!force){
                stop("Stopping analysis: ", directory, " already exists! Use force=TRUE to overwrite.")
            }
        } else {
            dir.create(directory, recursive=TRUE)
        }

        if(nrow(fc@module) > 0){
            write.table(fc@module, file.path(directory, "module.tsv"), sep="\t", row.names=FALSE)

            mean_summary <- mod_summary(fc, "mean")
            write.table(mean_summary, file.path(directory, "summary_mean.tsv"), sep="\t", row.names=FALSE)

            median_summary <- mod_summary(fc, "median")
            write.table(median_summary, file.path(directory, "summary_median.tsv"), sep="\t", row.names=FALSE)

            eg_summary <- mod_summary(fc, "eigengene")
            write.table(eg_summary, file.path(directory, "summary_eigengene.tsv"), sep="\t", row.names=FALSE)

            module_to_gmt(fc, directory=directory)
        }

        expr_f <- expr_data(fc, filter=fc@parameters$filter,
                            apply_vst=fc@parameters$apply_vst)
        selected <- select_genes(expr_f)
        writeLines(selected, file.path(directory, "selected_genes.txt"))

        if(length(fc@enrichment) > 0){
            for (stat in names(fc@enrichment)) {
                write.table(fc@enrichment[[stat]],
                            file.path(directory, paste0("enrichment_", stat, ".tsv")),
                            sep="\t", row.names=FALSE)
            }
        }

        if(nrow(fc@ora) > 0){
            write.table(fc@ora, file.path(directory, "ora.tsv"), sep="\t", row.names=FALSE)
        }

        if(length(fc@interactions) > 0){
            mod_ints <- lapply(names(fc@interactions), function(x){
                mod_int <- igraph::get.edgelist(fc@interactions[[x]])
                if(nrow(mod_int) > 0 ){
                    cbind(x, mod_int)
                }
            })
            int_df <- do.call("rbind", mod_ints)
            colnames(int_df) <- c("Module", "Gene1", "Gene2")
            write.table(int_df, file.path(directory, "interactions.tsv"), sep="\t", row.names=FALSE)
        }

        if(length(fc@parameters) > 0){
            params <- fc@parameters
            param_df <- data.frame(Parameter=names(params), Value=as.character(params))
            write.table(param_df, file.path(directory, "parameters.tsv"), sep="\t", row.names=FALSE)
        }
    }
)



#' # Run module overrepresentation analysis
#' fc <- mod_ora(fc, gmt)
#' # Check results
#' head(ora_data(fc))
#'
#' @rdname mod_ora
#' @export
setGeneric('mod_ora', function(fc, ...) {
  standardGeneric('mod_ora')
})

#' @rdname mod_ora
setMethod('mod_ora', signature('fcoex'),
          function(fc, gmt, verbose=FALSE) {
            #fc <- get_args(fc, vars=mget(ls()))
            if(!"gene" %in% names(gmt) | !"term" %in% names(gmt)){
              stop("The gmt object must contain two columns named 'term' and 'gene'")
            }
            if (verbose) {
              message('Running ORA')
              message("Using all genes in GMT file as universe.")
            }
            allgenes <- unique(gmt[, "gene"])
            if(is.null(fc@module_list)){
              warning("No modules in fcoex object! Did you run find_modules()?")
              return(fc)
            }
            mods <- fc@module_list
            res_list <- lapply(names(mods), ora, gmt, allgenes, mods)
            if (all(lapply(res_list, nrow) == 0)){
              warning("Enrichment is NULL. Either your gmt file is inadequate or your modules really aren't enriched for any of the pathways in the gmt file.")
              return(fc)
            }
            names(res_list) <- names(mods)
            
            res <- lapply(names(res_list), function(x){
              if(nrow(res_list[[x]]) > 0){
                as.data.frame(cbind(x, res_list[[x]]))
              }
            })
            res <- do.call(rbind, res)
            names(res)[names(res) == "x"] <- "Module"
            
            rownames(res) <- NULL
            fc@ora <- res
            return(fc)
          }
)

#' Retrieve over representation analysis (ORA) results
#'
#' @param fc Object of class \code{fcoex}
#'
#' @details This function returns the results of the \code{mod_ora} function on the
#' \code{fcoex} object. The ID column corresponds to pathways in the gmt file for which
#' genes in the modules were enriched. The Count column shows the number of genes in the
#' module that are enriched for each pathway. The GeneRatio column shows the proportion of
#' genes in the module enriched for a given pathway out of all the genes in the module
#' enriched for any given pathway. The BgRatio column shows the proportion of genes in a
#' given pathway out of all the genes in the gmt file. For more details, please refer to
#' the \code{clusterProfiler} package documentation.
#'
#' @return Object of class \code{data.frame} with ORA data
#'
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
          function(fc){
            return(fc@ora)
          })



#' Recluster cells based on fcoex module composition
#'
#' @param fc Object of class \code{fcoex}
#' @param hclust_method method for the hclust function. Defaults to "ward.D2".
#' @param dist_method  method for the dist function. Defaults to "manhattan".
#' @param k desired number of clustes. Defaults to 2.
#' @export
#' @rdname recluster
setGeneric("recluster", function(fc, ...) {
  standardGeneric("recluster")
})
#' @rdname recluster
setMethod("recluster", signature("fcoex"),
          function(fc,
                   hclust_method = "ward.D2",
                   dist_method = 'manhattan',
                   k = 2
){
            mod_idents <-list()
            for (i in names(fc@module_list)){
              print(i)
              expression_table <- fc@expression[fc@module_list[[i]],]
              d <- dist(t(as.matrix(expression_table)), method = dist_method)
              hc <- hclust(d, method = hclust_method)
              idents <- as.factor(cutree(hc, k))
              mod_idents[[i]] <- idents
          }
            fc@mod_idents <- mod_idents
            return(fc)
          })