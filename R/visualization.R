# Many parts of this code were directly adapted from the
# CEMiTool package. This include chunks of copied and pasted code
# located inside the source code for CEMiTool functions. 
# Functions that contained adapted code were explicit denoted.





#' @import ggplot2
#' @importFrom sna gplot.layout.fruchtermanreingold
#' @importFrom data.table melt
#' @importFrom ggrepel geom_label_repel
#' @importFrom igraph degree set_vertex_attr graph.adjacency
#' @import intergraph
#' @importFrom scales squish
#' @import stringr
#' @importFrom network as.matrix.network.adjacency as.matrix.network.edgelist get.vertex.attribute
#' @import grid
#' @importFrom pathwayPCA read_gmt
NULL



#' Network visualization
#'
#' Creates network visualizations based on the adjacency matrix
#' obtained with the find_cbf_modules method 
#'
#' This function was copied and adapted from the CEMiTool package.
#' The visualization of networks in this function is derivative of the 
#' intelectual work of CEMiTool's authors.
#' 
#' @param fc Object of class \code{fcoex}.
#' @param n number of nodes to label
#' @param min_elements Minimum number of elements in a module for it to be plotted. Defaults to 5.
#' @param ... Optional parameters.
#' @examples 
#' library(SingleCellExperiment) 
#' data("mini_pbmc3k")
#' targets <- colData(mini_pbmc3k)$clusters
#' exprs <- as.data.frame(assay(mini_pbmc3k, "logcounts"))
#' fc <- new_fcoex(exprs, targets)
#' fc <- discretize(fc)
#' fc <- find_cbf_modules(fc)
#' fc <- get_nets(fc)
#' @return Object of class \code{fcoex} with profile plots
#'
#' @rdname get_nets
#' @export
setGeneric('get_nets', function(fc,
                                         n = 10,
                                         min_elements = 5,
                                         ...) {
  standardGeneric('get_nets')
})

#' @rdname get_nets
setMethod('get_nets', signature('fcoex'),
          function(fc,
                   n = 10,
                   min_elements = 5,
                   ...) {
            if (length(fc@module_list) == 0) {
              stop("No modules in fcoex object! Did you run find_cbf_modules()?")
            }
            #fc <- get_args(fc, vars=mget(ls()))
            fc <- mod_colors(fc)
            module_cols <- fc@mod_colors
            mod_names <- names(fc@module_list)
            adjacency_full <- fc@adjacency_trimmed
            adj <- fc@adjacency_trimmed[, -1]
            rownames(adj) <- colnames(adj)
            res <- lapply(mod_names, function(name) {
              members_of_module <- fc@module_list[[name]]
              if (length(members_of_module) >= min_elements) {
                adj <- adj[members_of_module, members_of_module]
                adj <- as.matrix(adj)
                .plot_one_interaction(adj,
                                      n = n,
                                      color = module_cols[name],
                                      name = name)
              }
            })
            names(res) <- mod_names
            fc@coex_network_plot <- res[mod_names]
            return(fc)
          })

#' Set module colors mod_colors attribute
#' @param fc Object of class \code{fcoex}
#'
#' @return A vector with color names.
#' @examples 
#' data("fc")
#' mod_colors(fc)
#' @rdname mod_colors
#' @export
setGeneric("mod_colors", function(fc) {
  standardGeneric("mod_colors")
})

#' @rdname mod_colors
setMethod("mod_colors", signature("fcoex"),
          function(fc) {
            mod_names <- names(fc@module_list)
            nmod <- length(mod_names)
            cols <- fc@mod_colors
            if (nmod != 0) {
              if (length(fc@mod_colors) == 0) {
                if (nmod <= 16) {
                  cols <- rainbow(16, s = 1, v = 0.7)[seq_len(nmod)]
                } else {
                  cols <- rep(rainbow(16, s = 1, v = 0.7), ceiling(nmod / 16))[seq_len(nmod)]
                }
                names(cols) <- mod_names
              } else {
                if (is.null(names(fc@mod_colors))) {
                  warning("mod_colors should be a character vector with names corresponding to the modules")
                } else if (!all(sort(names(fc@mod_colors)) == sort(mod_names))) {
                  warning("mod_colors names do not match with modules!")
                }
              }
            }
            fc@mod_colors <- cols
            return(fc)
          })

#' Network visualization
#'
#' Creates a graph based on interactions provided
#' This function was copied and adapted from the CEMiTool package.
#' The visualization of networks in this function is derivative of the 
#' intelectual work of CEMiTool's authors.
#' 
#' 
#' @param adjacency_matrix An adajcency matrix from the \code{fcoex} object.
#' @param n Number of genes to be shown
#' @param color Color of the module to be plotted
#' @param name Name of the module to be plotted
#' @param ... Optional parameters.
#' @rdname plot_one_interaction
#' @return  A ggplot2 ('gg') object
.plot_one_interaction <- function(adjacency_matrix, n, color, name) {
  
  # comments below were also retained from the orignal CEMiTool code
  
  adj <- as.matrix(adjacency_matrix)
  ig_obj <- graph.adjacency(adj, weighted = TRUE, diag = FALSE)
  degrees <- igraph::degree(ig_obj, normalized = FALSE)
  ig_obj <-
    igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
  net_obj <- intergraph::asNetwork(ig_obj)
  m <-
    network::as.matrix.network.adjacency(net_obj) # get sociomatrix
  # get coordinates from Fruchterman and Reingold's force-directed plafcent algorithm.
  plotcord <-
    data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
  # or get it them from Kamada-Kawai's algorithm:
  # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  colnames(plotcord) <- c("X1", "X2")
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <-
    data.frame(plotcord[edglist[, 1],], plotcord[edglist[, 2],])
  plotcord$vertex.names <-
    as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  plotcord$Degree <-
    network::get.vertex.attribute(net_obj, "degree")
  plotcord[, "shouldLabel"] <- FALSE
  max_n <- min(n, length(degrees))
  int_hubs <- names(sort(degrees, decreasing = TRUE))[seq_len(max_n)]
  int_bool <- plotcord[, "vertex.names"] %in% int_hubs
  sel_vertex <- int_hubs
  colnames(edges) <-  c("X1", "Y1", "X2", "Y2")
  #edges$midX  <- (edges$X1 + edges$X2) / 2
  #edges$midY  <- (edges$Y1 + edges$Y2) / 2
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <-
    TRUE
  plotcord$Degree_cut <-
    cut(plotcord$Degree,
        breaks = 3,
        labels = FALSE)
  plotcord$in_mod <- TRUE
  pl <- ggplot(plotcord)  +
    geom_segment(
      data = edges,
      aes_(
        x =  ~ X1,
        y =  ~ Y1,
        xend =  ~ X2,
        yend =  ~ Y2
      ),
      size = 0.5,
      alpha = 0.8,
      colour = "#9e9e9e"
    ) +
    geom_point(aes_(
      x =  ~ X1,
      y =  ~ X2,
      size =  ~ Degree,
      alpha =  ~ Degree
    ), color = color) +
    geom_label_repel(
      aes_(
        x =  ~ X1,
        y =  ~ X2,
        label =  ~ vertex.names
      ),
      box.padding = unit(1, "lines"),
      data = function(x) {
        x[x$shouldLabel,]
      }
    ) +
    scale_colour_manual(values = c("#005E87")) +
    labs(title = name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white",
                                               colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  return(pl)
}


#' Retrieve fcoex net plots
#'
#' @param fc Object of class \code{fcoex}.
#' @return A plot corresponding to a fcoex analysis
#' @examples 
#' data("fc")
#' show_net(fc)
#' @rdname show_net
#' @export
setGeneric('show_net', function(fc) {
  standardGeneric('show_net')
})

#' @rdname show_net
setMethod('show_net', signature('fcoex'),
          function(fc) {
            return(fc@coex_network_plot)
            #       }
          })

#' ORA visualization
#'
#' Creates a bar plot with the results of module overrepresentation analysis
#' 
#' This function was copied and adapted from the CEMiTool package.
#' The visualization in this function is derivative of the 
#' intelectual work of CEMiTool's authors.
#' 
#'
#' @param fc Object of class \code{fcoex}.
#' @param n number of enrichments to show
#' @param pv_cut p-value significance cutoff. Default is 0.05.
#' @param ... parameters to plot_ora_single
#'
#' @return Object of class \code{fcoex} with ORA plots
#' @examples 
#' data("fc")
#' gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
#' gmt_in <- pathwayPCA::read_gmt(gmt_fname)
#' fc <- mod_ora(fc, gmt_in)
#' fc <- plot_ora(fc)
#' @rdname plot_ora
#' @export
setGeneric('plot_ora', function(fc, n = 10, pv_cut = 0.05, ...) {
  standardGeneric('plot_ora')
})

#' @rdname plot_ora
setMethod('plot_ora', signature('fcoex'),
          function(fc, n = 10, pv_cut = 0.05, ...) {
            if (length(fc@module_list) == 0) {
              stop("No modules in fcoex object! Did you run find_modules()?")
            }
            if (nrow(fc@ora) == 0) {
              stop("No ORA data! Did you run mod_ora()?")
            }
            
            #fc <- get_args(fc=fc, vars=mget(ls()))
            ora_splitted <- split(fc@ora, fc@ora$Module)
            module_cols <- fc@mod_colors
            res <- lapply(ora_splitted, function(x) {
              plot_ora_single(
                head(x, n = n),
                pv_cut = pv_cut,
                graph_color = module_cols[unique(x$Module)],
                title = unique(x$Module)
              )
            })
            fc@barplot_ora <- res
            return(fc)
          })


#' Retrieve fcoex ora plots
#'
#' This function was copied and adapted from the CEMiTool package.
#' The visualization in this function is derivative of the 
#' intelectual work of CEMiTool's authors.
#' 
#' @param fc Object of class \code{fcoex}.
#' @return A plot corresponding to a fcoex analysis
#' @examples 
#' data("fc")
#' show_ora(fc)
#' @rdname show_ora
#' @export
setGeneric('show_ora', function(fc) {
  standardGeneric('show_ora')
})

#' @rdname show_ora
setMethod('show_ora', signature('fcoex'),
          function(fc) {
            return(fc@barplot_ora)
            #       }
          })

#' ORA visualization for one module
#'
#'
#' This function was copied and adapted from the CEMiTool package.
#' The visualization in this function is derivative of the 
#' intelectual work of CEMiTool's authors.
#' 
#' @keywords internal
#'
#' @param es a data.frame from ora function containing only one module
#' @param ordr_by column to order the data.frame
#' @param max_length max length of a gene set name
#' @param pv_cut p-value cuttoff
#' @param graph_color color of bars
#' @param title title of the graph
#'
#' @return a list with ggplot2 object and the number of significant gene sets
plot_ora_single <-
  function(es,
           ordr_by = 'p.adjust',
           max_length = 50,
           pv_cut = 0.05,
           graph_color = "#4169E1",
           title = "Over Representation Analysis") {
    
    # comments below were also retained from the orignal CEMiTool code
    
    comsub <- function(x) {
      #split the first and last element by character
      d_x <- strsplit(x[c(1, length(x))], "")
      #search for the first not common element, and so, get the last matching one
      der_com <- match(FALSE, do.call("==", d_x)) - 1
      return(substr(x, 1, der_com + 1))
    }
    
    es[, "GeneSet"] <- es[, "ID"]
    
    # limits name length
    ovf_rows <- which(nchar(es[, "GeneSet"]) > max_length) # overflow
    ovf_data <- es[ovf_rows, "GeneSet"]
    test <- strtrim(ovf_data, max_length)
    dupes <- duplicated(test) | duplicated(test, fromLast = TRUE)
    if (sum(dupes) > 0) {
      test[dupes] <- ovf_data[dupes]
      test[dupes] <- comsub(test[dupes])
      max_length <- max(nchar(test))
    }
    
    es[ovf_rows, "GeneSet"] <-
      paste0(strtrim(test, max_length), "...")
    es[, "GeneSet"] <- stringr::str_wrap(es[, "GeneSet"], width = 20)
    
    # order bars
    lvls <- es[order(es[, ordr_by], decreasing = TRUE), "GeneSet"]
    es[, "GeneSet"] <- factor(es[, "GeneSet"], levels = lvls)
    
    es[, "alpha"] <- 1
    es[es[, ordr_by] > pv_cut, "alpha"] <- 0
    
    # Avoid 0's
    es[es[, ordr_by] > 0.8, ordr_by] <- 0.8
    my_squish <- function(...) {
      return(scales::squish(..., only.finite = FALSE))
    }
    
    # plot
    y_axis <- paste('-log10(', ordr_by, ')')
    pl <-
      ggplot(es,
             aes_string(
               x = "GeneSet",
               y = y_axis,
               alpha = "alpha",
               fill = y_axis
             )) +
      geom_bar(stat = "identity") +
      theme(axis.text = element_text(size = 8), legend.title = element_blank()) +
      coord_flip() +
      scale_alpha(range = c(0.4, 1), guide = "none") +
      labs(y = "-log10(adjusted p-value)", title = title, x = "") +
      geom_hline(
        yintercept = -log10(pv_cut),
        colour = "grey",
        linetype = "longdash"
      ) +
      scale_fill_gradient(
        low = "gray",
        high = graph_color,
        limits = c(2, 5),
        oob = my_squish
      )
    res <-
      list('pl' = pl, numsig = sum(es[, ordr_by] < pv_cut, na.rm = TRUE))
    return(res)
  }

#' @title
#' Save fcoex object plots
#'
#' @description
#' Save plots into the directory specified by the \code{directory} argument.
#' Note: If no directory is specified, it will save to tempdir().
#' A possible option is setting directory = "./Plots"
#' 
#' This function was modified from the CEMiTool package.
#' Chunks of code were retained "as is"
#'
#' @param fc Object of class \code{fcoex}.
#' @param name The name of the file to be saved.
#' @param directory Directory into which the files will be saved.
#' @param force If the directory exists, execution will not stop.
#' @return A pdf file or files with the desired plot(s)
#' @examples 
#' data(fc)
#' save_plots(fc, name = "Example")
#' @rdname save_plots
#' @export
setGeneric('save_plots', function(fc, name,
                                  force = FALSE,
                                  directory = "tempdir()") {
  standardGeneric('save_plots')
})

#' @rdname save_plots
setMethod('save_plots', signature('fcoex'),
          function(fc,
                   name,
                   force = FALSE,
                   directory = "tempdir()") {
            if (dir.exists(directory)) {
              if (!force) {
                stop("Stopping analysis: ",
                     directory,
                     " already exists! Use force=TRUE to overwrite.")
              }
            } else{
              dir.create(directory, recursive = TRUE)
            }
            plots <- list(fc@coex_network_plot, fc@barplot_ora)
            all_plots <- c("interaction", "ora")
            names(plots) <- all_plots
            plots <- Filter(function(x)
              length(x) >= 1, plots)
            if (length(plots) < 2) {
              message(
                "Saving available plots."
              )
            }
            lapply(names(plots), function(pl) {
              pdf(file = file.path(directory, paste0(name, "_", pl, ".pdf")))
              print(plots[[pl]])
              dev.off()
            })
            
          })
