#' removeEvidenceFromPed: a function for removing evidence from specific individuals in a ped object.
#'
#' @param pped A ped object with information of the genotyped members. The ped object must be in Familias format.
#' @param idNotEv A set of individuals whom evidence should be removed. 
#' @import paramlink
#' @import igraph
#' @export
#' @return A ped object.

removeEvidenceFromPed<-function(pped,idNotEv){
 for(i in 1:pped$nMark){
  pped$markerdata[[i]][idNotEv,]<-c(0,0)
 }
 pped$available <- pped$orig.ids[!pped$orig.ids%in%idNotEv]
 return(pped)
}
