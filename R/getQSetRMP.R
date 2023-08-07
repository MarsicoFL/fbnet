#' getGenotypeTables: a function for obtaining genotypetables after variable elimination and using available genetic evidence.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param lqp list of individuals genotypes.
#' @import paramlink
#' @import igraph
#' @export
#' @return A dataframe with genotype probabilities.

getQSetRMP<-function(bn,lqp){
 count<-rmp<-1
 for(i in seq_along(lqp)){
     sys <- unlist(lapply(strsplit(names(lqp)[i],"_"),function(x){return(x[2])}))
     rmp <- rmp * prod(bn$alelFreq[[sys]][lqp[[i]]]) * 2^(lqp[[i]][1]!=lqp[[i]][2])
 }
 return(rmp)
}
