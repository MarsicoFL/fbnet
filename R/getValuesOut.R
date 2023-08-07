#' getValuesOut: a function for getting out variables with zero probability in the bayesian network
#'
#' @param cpt conditional probability table from the bayesian network
#' @param condVar variables from the conditioning table
#' @import paramlink
#' @import igraph
#' @export
#' @return A processed conditional probability table 

getValuesOut<-function(cpt,condVar=c()){
   ccol<-colnames(cpt)[!colnames(cpt)%in%c(condVar,"prob")]
   loutvarcpt<-list()
   for(j in seq_along(ccol)){
       cpt0<-unique(cpt[cpt[,"prob"]==0,ccol[j]])
       cpt1<-unique(cpt[cpt[,"prob"]!=0,ccol[j]])
       valuesOut<-setdiff(cpt0,cpt1)
       if(length(valuesOut)>0){
        voutName <- ccol[j]
        loutvarcpt[[voutName]] <-c(loutvarcpt[[voutName]],valuesOut)
       }
   }
   return(loutvarcpt)
}

