#' setOrdering: a function for selecting the ordering method in the elimination process.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param ordMethod Ordering method.
#' @param vars Vars
#' @param orderElim Order elimination criteria.
#' @import paramlink
#' @import igraph
#' @export
#' @return A bayesian network after ordering process.

setOrdering<-function(bn,ordMethod,vars=NULL,orderElim=NULL){

 bn$ordMethod <- ordMethod

 if(ordMethod=="id"){
  if(!is.null(vars)){
    bnordering<-names(bn$vars)[!names(bn$vars)%in%bn$Q]
    bnordering<-bn[["ordering"]][bn[["ordering"]]%in%vars]
  }else{
   bnordering<-names(bn$vars)[!names(bn$vars)%in%bn$Q]
  }
 }

 if(ordMethod%in%c("min_degree","min_fill")){
   bnordering<-minOrdering(bn,vars,ordMethod)
 }

 if(ordMethod=="fixed"){
  if(is.null(orderElim)) warning("ordMethod=fixed, but no orderElim provided.\n")
  if(!all(orderElim%in%names(bn$CPTs))) warning("check orderElim bariable names.\n")
  bnordering<-orderElim
 }

 return(bnordering)
}

