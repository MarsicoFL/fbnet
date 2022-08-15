#' evidencePrunning: a fuction for pruning instantiated variables.
#'
#' @param bn A bayesian network (output of buildBN function).
#' @import paramlink
#' @import igraph
#' @export
#' @return A preprocessed bayesian network.

evidencePrunning<-function(bn){
   E<-bn$E
   for(i in seq_along(E)){
     evidencia<-names(E)[i]
     nei <- neighbors(bn$DAG,evidencia,mode="out")
     if(length(nei)>0){
      for(j in seq_along(nei)){
       target   <-names(nei)[j]
       cpt <- bn$CPTs[[target]]


       cpt <- cpt[ cpt[,evidencia]%in%E[[i]],]

       bn$CPTs[[target]]<-cpt

      }
     }
     cpt <- bn$CPTs[[ evidencia ]]
     cpt <- cpt[ cpt[,evidencia]%in%E[[i]],]

     bn$CPTs[[evidencia]] <- cpt
    }

    return(bn)

  }
