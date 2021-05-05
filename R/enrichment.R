# Many parts of this code were directly adapted from the
# CEMiTool package. This include chunks of copied and pasted code
# located inside the source code for CEMiTool functions.
# Functions that contained adapted code were explicit denoted.

# Imports --------

#' @importFrom clusterProfiler enricher
NULL

# Set Over Representation Analysis (ORA) --------

#'  Over Representation Analysis (ORA)
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
setGeneric("mod_ora", function(fc, gmt, verbose = FALSE) {
  standardGeneric("mod_ora")
})

#' @rdname mod_ora
setMethod(
  "mod_ora", signature("fcoex"),
  function(fc, gmt, verbose = FALSE) {
    if (!is(gmt, "pathwayCollection")) {
      stop(
        "The gmt object must be loaded via pathwayPCA::read_gmt and be a pathwayCollectionobject"
      )
    }
    
    pathwayPCA_gmt <- gmt
    gmt_df <- pathwayPCA_gmt$pathways
    
    gmt_df <-
      as.data.frame(unlist(gmt_df), use.names = TRUE)
    pathway_sizes <-
      unlist(lapply(pathwayPCA_gmt$pathways, length))
    colnames(gmt_df) <- "gene"
    gmt_df$term <- rep(pathwayPCA_gmt$TERMS, pathway_sizes)
    
    
    gmt_df <- gmt_df[, c("term", "gene")]
    
    
    if (verbose) {
      message("Running ORA")
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
  }
)


# Performs Over Representation Analysis for a list of genes and a GMT
# This function was modified from the CEMiTool package.
# Chunks of code were retained "as is"
# @keywords internal
#
# @param topgenes a vector of genes
# @param gmt.list a gmt from prepare.gmt function
# @param allgenes a vector containing all genes to be considered as universe
#
# @return a data.frame containing the results
#
#
ora <- function(mod_name, gmt_list, allgenes, mods) {
  if (missing(allgenes)) {
    message("Using all genes in GMT file as universe.")
    allgenes <- unique(gmt_list[, "gene"])
  }
  topgenes <- mods[[mod_name]]
  enriched <- clusterProfiler::enricher(
    gene = topgenes,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    universe = allgenes,
    TERM2GENE = gmt_list
  )
  # TERM2NAME = gmt_list[['term2name']])
  if (!is.null(enriched) && !is.logical(enriched)) {
    result <- enriched@result
  } else {
    if (mod_name != "Not.Correlated") {
      warning("Enrichment for module ", mod_name, " is NULL")
    }
    result <- data.frame(
      Module = character(),
      ID = character(),
      Description = character(),
      GeneRatio = numeric(),
      BgRatio = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = numeric(),
      stringsAsFactors = FALSE
    )
  }
  return(result)
}
