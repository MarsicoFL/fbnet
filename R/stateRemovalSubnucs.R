#' stateRemovalSubnucs: a fuctiong for variable state pruning.
#'
#' @param bn A bayesian network (output of buildBN function).
#' @param verbose Computation output.
#' @import paramlink
#' @import igraph
#' @export
#' @return A preprocessed bayesian network.

stateRemovalSubnucs<-function(bn,verbose=FALSE){
cpt1 <- cpt2 <- NULL

  for(i in seq_along(bn$ped$subnucs)){
  if(bn$ped$subnucs$offspring){
   bby<-intersect(colnames(cpt1),colnames(cpt2))
   bby<-bby[!bby%in%"prob"]
   res<-merge(cpt1,cpt2,by=bby)
   iprob<-grep("prob",colnames(res))
   res<-cbind(res[,-iprob],prob=apply(res[,iprob],1,prod))
  }
 }
}

