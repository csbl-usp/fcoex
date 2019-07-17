#' @importFrom clusterProfiler enricher
NULL


#' Read a GMT file
#'
#' copied ipsis litteris from CEMiTool package.
#' @param fname GMT file name.
#' @return A list containing genes and description of each pathway
#' @examples
#' # Read example gmt file
#' gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
#' gmt_in <- read_gmt(gmt_fname)
#'
#' @export

read_gmt <- function(fname){
  res <- list(genes=list(), desc=list())
  gmt <- file(fname)
  gmt_lines <- readLines(gmt)
  close(gmt)
  gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, split="\t")))
  gmt_names <- sapply(gmt_list, '[', 1)
  gmt_desc <- lapply(gmt_list, '[', 2)
  gmt_genes <- lapply(gmt_list, function(x){x[3:length(x)]})
  names(gmt_desc) <- names(gmt_genes) <- gmt_names
  res <- do.call(rbind, lapply(names(gmt_genes),
                               function(n) cbind.data.frame(term=n, gene=gmt_genes[[n]], stringsAsFactors=FALSE)))
  res$term <- as.factor(res$term)
  return(res)
}

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
ora <- function(mod_name, gmt_list, allgenes, mods){
  if(missing(allgenes)) {
    message("Using all genes in GMT file as universe.")
    allgenes <- unique(gmt_list[, "gene"])
  }
  topgenes <- mods[[mod_name]]
  enriched <- clusterProfiler::enricher(gene = topgenes,
                                        pvalueCutoff = 1,
                                        qvalueCutoff = 1,
                                        universe = allgenes,
                                        TERM2GENE = gmt_list)
  # TERM2NAME = gmt_list[['term2name']])
  if (!is.null(enriched) && !is.logical(enriched)) {
    result <- enriched@result
  } else {
    if(mod_name != "Not.Correlated"){
      warning("Enrichment for module ", mod_name, " is NULL")
    }
    result <- data.frame(Module=character(), ID=character(),
                         Description=character(),
                         GeneRatio=numeric(), BgRatio=numeric(),
                         pvalue=numeric(), p.adjust=numeric(),
                         qvalue=numeric(), geneID=character(),
                         Count=numeric(), stringsAsFactors=FALSE)
  }
  return(result)
}

