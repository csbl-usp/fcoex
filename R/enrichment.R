#' @importFrom clusterProfiler enricher
NULL

# Performs Over Representation Analysis for a list of genes and a GMT
#
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

  pathwayPCA_gmt <- gmt_list
  gmt_df <- pathwayPCA_gmt$pathways
  names(gmt_df) <- pathwayPCA_gmt$TERMS
  
  gmt_df <- as.data.frame(unlist(gmt_df))
  colnames(gmt_df) <- "gene"
  gmt_df$term <- rownames(gmt_df)
  
  
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
