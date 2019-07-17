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
NULL


#' Network visualization
#'
#' Creates a graph based on interactions provided
#'
#' @param fc Object of class \code{fcoex}.
#' @param n number of nodes to label
#' @param min_elements Minimum number of elements in a module for it to be plotted. Defaults to 5. 
#' @param ... Optional parameters.
#'
#' @return Object of class \code{fcoex} with profile plots
#'
#' @rdname plot_interactions
#' @export
setGeneric('plot_interactions', function(fc, n=10, min_elements = 5,...) {
    standardGeneric('plot_interactions')
})

#' @rdname plot_interactions
setMethod('plot_interactions', signature('fcoex'),
    function(fc, n=10, min_elements = 5, ...) {
        if(length(fc@module_list) == 0){
            stop("No modules in fcoex object! Did you run find_cbf_modules()?")
        }
        #fc <- get_args(fc, vars=mget(ls()))
        fc <- mod_colors(fc)
        module_cols <- fc@mod_colors
        mod_names <- names(fc@module_list)
        adjacency_full <- fc@adjacency
        adj <- fc@adjacency[,-1]
        rownames(adj) <- colnames(adj)
        res <- lapply(mod_names, function(name){
          members_of_module <- fc@module_list[[name]]
                  if(length(members_of_module) >= min_elements){
                  adj <- adj[members_of_module,members_of_module]
                  adj <- as.matrix(adj)
                  .plot_one_interaction(adj,
                                    n=n, color=module_cols[name], name=name)
                  }
               })
        names(res) <- mod_names
        fc@interaction_plot <- res[mod_names]
        return(fc)
    })


#' Network visualization
#'
#' Creates a graph based on interactions provided
#'
#' @param fc Object of class \code{fcoex}.
#' @param ... Optional parameters.
.plot_one_interaction <- function(adjacency_matrix, n, color, name){
    adj <- as.matrix(adjacency_matrix)
    ig_obj <- graph.adjacency(adj, weighted = T, diag=F)
    degrees <- igraph::degree(ig_obj, normalized = FALSE)
    ig_obj <-
      igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
    net_obj <- intergraph::asNetwork(ig_obj)
    m <-      network::as.matrix.network.adjacency(net_obj) # get sociomatrix
    # get coordinates from Fruchterman and Reingold's force-directed plafcent algorithm.
    plotcord <-
      data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
    # or get it them from Kamada-Kawai's algorithm:
    # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
    colnames(plotcord) <- c("X1", "X2")
    edglist <- network::as.matrix.network.edgelist(net_obj)
    edges <-
      data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
    plotcord$vertex.names <-
      as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
    plotcord$Degree <-
      network::get.vertex.attribute(net_obj, "degree")
    plotcord[, "shouldLabel"] <- FALSE
    max_n <- min(n, length(degrees))
    int_hubs <- names(sort(degrees, decreasing = TRUE))[1:max_n]
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
        geom_segment(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2),
                     size = 0.5, alpha=0.8, colour="#9e9e9e") +
        geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree), color=color) +
      geom_label_repel(
        aes_( x =  ~ X1,   y =  ~ X2,  label =  ~ vertex.names        ),
                         box.padding=unit(1, "lines"),
                         data=function(x){x[x$shouldLabel, ]}) +
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
        panel.grid = ggplot2::element_blank())

    return(pl)
}


#' Retrieve fcoex object plots
#'
#' @param fc Object of class \code{fcoex}.
#' @return A plot corresponding to a fcoex analysis
#' @rdname show_plot
#' @export
setGeneric('show_plot', function(fc, value) {
    standardGeneric('show_plot')
})

#' @rdname show_plot
setMethod('show_plot', signature('fcoex'),
    function(fc 
             ) {
            return(fc@interaction_plot)
#       }
    })

#' ORA visualization
#'
#' Creates a bar plot with the results of module overrepresentation analysis
#'
#' @param fc Object of class \code{fcoex}.
#' @param n number of modules to show
#' @param pv_cut p-value significance cutoff. Default is 0.05.
#' @param ... parameters to plot_ora_single
#'
#' @return Object of class \code{fcoex} with ORA plots
#' @rdname plot_ora
#' @export
setGeneric('plot_ora', function(fc, ...) {
  standardGeneric('plot_ora')
})

#' @rdname plot_ora
setMethod('plot_ora', signature('fcoex'),
          function(fc, n=10, pv_cut=0.05, ...){
            if(length(fc@module_list) == 0){
              stop("No modules in fcoex object! Did you run find_modules()?")
            }
            if(nrow(fc@ora) == 0){
              stop("No ORA data! Did you run mod_ora()?")
            }
            
            #fc <- get_args(fc=fc, vars=mget(ls()))
            ora_splitted <- split(fc@ora, fc@ora$Module)
            module_cols <- fc@mod_colors
            res <- lapply(ora_splitted, function(x){
              plot_ora_single(head(x, n=n),
                              pv_cut=pv_cut,
                              graph_color=module_cols[unique(x$Module)],
                              title=unique(x$Module),
                              ...)
            })
            modules <- names(res)
            modules <- modules[order(as.numeric(stringr::str_extract(modules, "\\d+")))]
            fc@barplot_ora <- res[modules]
            return(fc)
          })

#' ORA visualization for one module
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
plot_ora_single <- function(es, ordr_by='p.adjust', max_length=50, pv_cut=0.05,
                            graph_color="#4169E1", title="Over Representation Analysis"){
  
  comsub <- function(x){
    #split the first and last element by character
    d_x <- strsplit(x[c(1, length(x))], "")
    #search for the first not common element, and so, get the last matching one
    der_com <- match(FALSE, do.call("==", d_x))-1
    return(substr(x, 1, der_com + 1))
  }
  
  es[, "GeneSet"] <- es[, "ID"]
  
  # limits name length
  ovf_rows <- which(nchar(es[, "GeneSet"]) > max_length) # overflow
  ovf_data <- es[ovf_rows, "GeneSet"]
  test <- strtrim(ovf_data, max_length)
  dupes <- duplicated(test) | duplicated(test, fromLast=TRUE)
  if(sum(dupes) > 0){
    test[dupes] <- ovf_data[dupes]
    test[dupes] <- comsub(test[dupes])
    max_length <- max(nchar(test))
  }
  
  es[ovf_rows, "GeneSet"] <-  paste0(strtrim(test, max_length), "...")
  es[, "GeneSet"] <- stringr::str_wrap(es[, "GeneSet"], width = 20)
  
  # order bars
  lvls <- es[order(es[, ordr_by], decreasing=TRUE), "GeneSet"]
  es[, "GeneSet"] <- factor(es[, "GeneSet"], levels=lvls)
  
  es[, "alpha"] <- 1
  es[es[, ordr_by] > pv_cut, "alpha"] <- 0
  
  # Avoid 0's
  es[es[, ordr_by] > 0.8, ordr_by] <- 0.8
  my_squish <- function(...){
    return(scales::squish(..., only.finite=FALSE))
  }
  
  # plot
  y_axis <- paste('-log10(', ordr_by, ')')
  pl <- ggplot(es, aes_string(x="GeneSet", y=y_axis, alpha="alpha", fill=y_axis)) +
    geom_bar(stat="identity") +
    theme(axis.text=element_text(size=8), legend.title=element_blank()) +
    coord_flip() +
    scale_alpha(range=c(0.4, 1), guide="none") +
    labs(y="-log10(adjusted p-value)", title=title, x="") +
    geom_hline(yintercept=-log10(pv_cut), colour="grey", linetype="longdash") +
    scale_fill_gradient(low="gray", high=graph_color, limits=c(2, 5), oob=my_squish)
  res <- list('pl'=pl, numsig=sum(es[, ordr_by] < pv_cut, na.rm=TRUE))
  return(res)
}

#' @title
#' Save fcoex object plots
#'
#' @description
#' Save plots into the directory specified by the \code{directory} argument.
#'
#' @param fc Object of class \code{fcoex}.
#' @param name The name of the file to be saved.
#' @param directory Directory into which the files will be saved.
#' @param force If the directory exists, execution will not stop.
#' @return A pdf file or files with the desired plot(s)
#' @rdname save_plots
#' @export
setGeneric('save_plots', function(fc, ...) {
  standardGeneric('save_plots')
})

#' @rdname save_plots
setMethod('save_plots', signature('fcoex'),
          function(fc, name, force=FALSE, directory="./Plots") {
            if(dir.exists(directory)){
              if(!force){
                stop("Stopping analysis: ", directory, " already exists! Use force=TRUE to overwrite.")
              }
            }else{
              dir.create(directory, recursive=TRUE)
            }
            plots <- list(fc@interaction_plot, fc@barplot_ora)
            all_plots<- c("interaction", "ora")
            names(plots) <- all_plots
            plots <- Filter(function(x) length(x) >= 1, plots)
            if(length(plots) < 2){
              message("Some plots have not been defined. Please run the appropriate plot functions. Saving available plots.")
            }
            lapply(names(plots), function(pl){
              pdf(file=file.path(directory, paste0(name,"_",pl,".pdf")))
              print(plots[[pl]])
              dev.off()
            })
            
          })
