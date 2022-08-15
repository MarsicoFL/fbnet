#' factorHeteroFounders: a function for multiplying probabilities in case of heterocigote founders.
#'
#' @param topn Format parameter.
#' @param resQ List of CPTs.
#' @import paramlink
#' @import igraph
#' @export
#' @return A MAP from the probability table.

getMAP<-function(resQ,topn=3){
  
 resQtop <- lapply(resQ,function(x){
                        x<-x[order(x[,"prob"],decreasing=TRUE)[1:min(topn,nrow(x))],]
                        return(x)
                    })}
