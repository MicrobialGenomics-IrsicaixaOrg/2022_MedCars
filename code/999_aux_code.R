#' From : https://github.com/rusher321/microbiotaPair/blob/master/R/network.R
#' sparcc network
#'
#' @param dat data.frame or matrix, row is sample ID
#' @param cutoff numberic, the sparcc coefficient of species-species indicate
#' the cutoff of eage
#' @param main character, the network title
#' @param count logistic, indicated the microbiota data is count of relative abundance
#' @param layout character, the paremeter from the igraph indicated the network layout
#' @param ...
#'
#' @return figure
#' @export
#'
#' @usage sparccNet(dat = microbiota, cutoff = 0.3, main = "test", count = F)
#' @examples
#'
#' library(dplyr)
#' data("physeq_data")
#' physeq <- physeq_data
#' microbitota <- otu_table(phyloseq)
#' microbitoFilter <- filterPer(microbiota, 2, 0.5)
#' sparccNet(dat = microbiotaFilter, cutoff = 0.3, main = "test", count = "F)
sparccNet <- function(dat, cutoff, main, count ,layout = "layout.circle",...){

  # library(igraph)
  # library(Matrix)
  # library(SpiecEasi)

  #
  renorm <- function(dat){
    # normlization function
    trans <- function(x){
      y <- x/(sum(x))
      return(y)
    }
    # tran the dataframe
    dat2 <- t(dat)
    dat3 <- t(apply(dat2, 2, trans))

    return(dat3)
  }

  # the function is count data or relative abundan
  if(count){
    datRenorm <- dat
  }else{
    datRenorm <- round(renorm(dat)*10^7)
  }
  basesparcc <- SpiecEasi::sparcc(datRenorm)
  graph <- basesparcc$Cor
  num.row <- nrow(graph)

  # tran the corelation to adj matrox

  for(i in 1:num.row){
    for(j in 1:num.row){
      a <- graph[i, j]
      graph[i,j] <- ifelse(abs(a) >= cutoff, a, 0)
    }
  }

  diag(graph) <- 0
  igraph <- SpiecEasi::adj2igraph(Matrix::Matrix(graph, sparse=TRUE))

  # set edge color，postive correlation
  # postive correlation is red, negative correlation is blue
  igraph.weight = E(igraph)$weight
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#fc9272", ifelse(E.color<0, "#31a354","grey"))
  E(igraph)$color = as.character(E.color)

  # add the edge width
  E(igraph)$width = abs(igraph.weight)*4
  # add the node color
  V(igraph)$color <- "#636363"
  plot(igraph,
       layout=layout,
       vertex.label=colnames(dat),
       main = main ,...)

}

