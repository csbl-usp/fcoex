
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
new_fcoex <-
  function(expression = data.frame(),
           target = vector()) {
    fc <- new("fcoex", expression = expression, target = target)
    msg <- "Created new fcoex object."
    message(msg)
    return(fc)
  }