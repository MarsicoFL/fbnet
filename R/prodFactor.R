#' prodFactor: a function for performing product between probability tables.
#'
#' @param laux probability distribution aux
#' @import paramlink
#' @import igraph
#' @export
#' @return A dataframe with probabilities.

prodFactor<-function(laux){

T1 <- expand.grid(list(A=1:4,B=1:2,C=1:3))
T2 <- expand.grid(list(B=1:5,C=1:2,D=1:3))
T3 <- expand.grid(list(A=c(2,3,5),B=1:3,C=1:3,E=1:2))
T1<-cbind(T1,prob=stats::runif(nrow(T1)))
T2<-cbind(T2,prob=stats::runif(nrow(T2)))
T3<-cbind(T3,prob=stats::runif(nrow(T3)))
laux<-list(T1=T1,T2=T2,T3=T3)


 if(length(laux)==1) return(laux[[1]])
 pprod <- laux[[1]]

 if(length(laux)>1){
  for(i in 2:length(laux)){
   ccol <- intersect(colnames(pprod),colnames(laux[[i]]))
   ccol <- ccol[!ccol%in%"prob"]
   pprod <- suppressWarnings(merge(pprod,laux[[i]],by=ccol))
  }
  iprob <- grep("prob",colnames(pprod))
  df <-cbind(pprod[,-iprob],prob=apply(pprod[,iprob],1,prod))
 }else{
  df <- pprod
 }

 return(df)
}

