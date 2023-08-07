#' minOrdering: a function for getting an ordering of bayesian network variables not in Q using min fill criteria on interaction graphs.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param method Elimination method, min_degree or min_fill
#' @param vars Subset of tables where the order is calculated
#' @import paramlink
#' @import igraph
#' @export
#' @return A bayesian network after ordering process.

minOrdering<-function(bn,vars=NULL,method=c("min_degree","min_fill")[1]){
 g <- make_empty_graph(directed=FALSE)

 if(is.null(vars)){
  lcpt <- bn$CPTs
 }else{
  lcpt <- bn$CPTs[vars]
 }

 for(i in seq_along(lcpt)){
      vcond<-colnames(lcpt[[i]])
      if(!is.null(vcond)){
       vcond<-vcond[!vcond%in%c("prob",bn$Q)]
       if(length(vcond)>0){
         gg <- make_full_graph(length(vcond))
         V(gg)$name <- vcond
         g <- g + gg
       }
      }
  }
  allvars <- V(g)$name

  #layout(matrix(1:5,1,5))
  #plot(g)
  ord <- c()
  for(i in 1:(vcount(g)-1)){
    if(method=="fill"){
     tc <- transitivity(g,"local") #porbabilidad que un par de vecinos esten conectados
     tc[is.na(tc)] <- 1000         #mando los que estan aislados (!) al final
     dg <- degree(g)
     numVecNoConectados<-(1-tc)*dg*(dg-1)/2
     imin <- which.min(numVecNoConectados)
    }else{
     imin<-which.min(degree(g))
    }
    ord <- c(ord,names(imin))
    nn  <- neighbors(g,imin)
    if(length(nn)>1){  #enlazo vecinos del que voy a eliminar
      gg <- make_full_graph(length(nn))
      V(gg)$name <- names(nn)
      g <- g + gg
    }
    g <- delete_vertices(g,names(imin))
    #plot(g)
  }
  ord <- c(ord,setdiff(allvars,ord))

  return(ord)
}

