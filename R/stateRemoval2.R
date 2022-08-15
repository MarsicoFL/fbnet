#' stateRemoval2: a function for processing the bayesian network. It implements another approach from the described in stateRemoval function.
#'
#' @param bn A bayesian network (output of buildBN function).
#' @param verbose Computation output.
#' @import paramlink
#' @import igraph
#' @export
#' @return A preprocessed bayesian network.

stateRemoval2 <- function(bn,verbose=FALSE){
 lCPTs<-bn$CPTs
 reverseSplit<-function (inList) {
    if (length(inList) == 0) {
        return(inList)
    }
    lens = sapply(inList, length)
    nms = rep(names(inList), lens)
    vals = unlist(inList)
    split(nms, vals)
}
 lWhichTable <- reverseSplit(lapply(lCPTs,function(x){colnames(x)[!colnames(x)%in%"prob"]}))


 checkList <- names(lCPTs)
 while(length(checkList)>0){
  i<-which(names(lCPTs)==checkList[1])
  cpt <- lCPTs[[i]]
  lout<-getValuesOut(cpt)
  if(length(lout)==0){
   checkList<-checkList[-1]
   if(verbose)cat("no",names(lCPTs)[i],length(checkList),"\n")
  }else{
   for(ilout in seq_along(lout)){
    variable <- names(lout)[ilout]
    #recorro todas las tablas sacando los valores de lout
    for(sTable in lWhichTable[[variable]]){
     ikeep<-which(!lCPTs[[sTable]][,variable]%in%lout[[ilout]])
     lCPTs[[sTable]] <- lCPTs[[sTable]][ikeep,]
     cnames <- colnames(lCPTs[[sTable]])
     checkList <- union(checkList,cnames[!cnames%in%c("prob")])
    }
   }
   if(verbose)cat("si",names(lCPTs)[i],length(checkList),"\n")
  }
 }
 bn$CPTs <- lCPTs
 return(bn)
}

