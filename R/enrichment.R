#' @importFrom data.table fread setDFs
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


#' Module Overrepresentation Analysis
#'
#' Performs overrepresentation analysis for each co-expression module found.
#'
#' @param fc Object of class \code{fcoex}.
#' @param gmt Object of class \code{data.frame} with 2 columns, one with
#' pathways and one with genes
#' @param verbose logical. Report analysis steps.
#' @param ... Optional parameters.
#'
#' @return Object of class \code{fcoex}
#'
#' @seealso \code{\link{ora_data}}
#'
#' @examples
#' # Get example fcoex object
#' data(fc)
#' # Read gmt file
 gmt <- read_gmt(system.file('extdata', 'pathways.gmt',
                 package='CEMiTool'))
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