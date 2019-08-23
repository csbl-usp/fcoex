#' Processed subset of the pbmc3k dataset from PBMC genomics
#'
#' A subset with 600 sampled cells and the top 1700 variable 
#' genes from the  TENxPBMCData package pbmc3k dataset. 
#' 
#' Preprocessed in accordance to OSCA (August 2019, https://osca.bioconductor.org/)
#' 
#' scater::normalized .
#' PCA with 50 components.
#' snn graph on the PCA space + louvain clustering to yield 8 clusters . 
#' UMAP already ran 
#'
#' @name mini_pbmc3k 
#' @docType data
#' @usage data(mini_pbmc3k)
#' @format An object of class \code{SingleCellExperiment}
#' @keywords datasets
#' @source \href{https://bioconductor.org/packages/release/data/experiment/html/TENxPBMCData.html}
#' @examples
#' data(mini_pbmc3k)
#' mini_pbmc3k
"mini_pbmc3k"


#' Example fcoex object
#'
#' Example fcoex object processed from the mini_pbmc3k dataset.
#' 
#' @name fc 
#' @docType data
#' @usage data(fc)
#' @format An object of class \code{fcoex}
#' @keywords datasets
#' @examples
#' data(fc)
#' fc
"fc"