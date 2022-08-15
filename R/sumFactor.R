#' prodFactor: a function for performing sum between probability tables.
#'
#' @param cpt Conditional probability table
#' @param Z factor
#' @import paramlink
#' @import igraph
#' @export
#' @return A dataframe with probabilities.

sumFactor<-function(cpt,Z){
 X <- colnames(cpt)
 Y <- X[!X%in%c(Z,"prob")]
 if(length(Y)==0) return(1)  #Z=Y!!
 

 if(length(Y)>1){
  f<-stats::aggregate(cpt[,"prob"],by=as.list(cpt[,Y]),sum)
 }else{ 
  f<-stats::aggregate(cpt[,"prob"],by=list(cpt[,Y]),sum)
 } 

 colnames(f)<-c(Y,"prob")
 return(f)
}
