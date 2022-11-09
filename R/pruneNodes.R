#' pruneNodes: a fuction for clasical pruning in bayesian networks.
#'
#' @param bn A bayesian network (output of buildBN function).
#' @import igraph
#' @export
#' @return A preprocessed bayesian network.

pruneNodes<-function(bn){
  E <- bn$E
  Q <- bn$Q  
  U <- union(paste0(rep(names(E),each=2),c("_p","_m"),sep=""),Q)
  dout    <- degree(bn$DAG,mode="out")
  leavesfb<-names(dout)[dout==0]
  out<-leavesfb[!leavesfb%in%U]
  while(length(out)>0){
   bn$DAG<-delete_vertices(bn$DAG,out)
   bn$CPTs[out] <- NULL
   
   dout    <- degree(bn$DAG,mode="out")
   leavesfb<-names(dout)[dout==0]
   
   out<-leavesfb[!leavesfb%in%U]
  }
  return(bn)
}

