#' stateRemoval: a function for processing the bayesian network.
#'
#' @param bn A bayesian network (output of buildBN function).
#' @import paramlink
#' @import igraph
#' @export
#' @return A preprocessed bayesian network.

stateRemoval <- function(bn){

   lCPT<-bn$CPTs
   iout <- grep("_pm",names(lCPT))
   lout<-list()
   for(i in seq_along(lCPT)){
     if(i%in%iout) next
     cole <-colnames(lCPT[[i]])[colnames(lCPT[[i]])%in%names(bn$E)]
     if(length(cole)>0 & ncol(lCPT[[i]])>2){
      cpt <- lCPT[[i]][apply(lCPT[[i]][,cole,drop=FALSE],1,function(x){all(x==bn$E[cole])}),]
      if(nrow(cpt)==0) warning("Incompatible data in ",names(lCPT)[i]," table\n")
      lCPT[[i]] <- cpt

      ccol<-colnames(cpt)[!colnames(cpt)%in%c(cole,"prob")]
      for(j in seq_along(ccol)){
       cpt0<-unique(cpt[cpt[,"prob"]==0,ccol[j]])
       cpt1<-unique(cpt[cpt[,"prob"]!=0,ccol[j]])
       valuesOut<-setdiff(cpt0,cpt1)
       if(length(valuesOut)>0){
        voutName <- ccol[j]
        lout[[voutName]] <-c(lout[[voutName]],valuesOut)
       }
      }
     }
   }

   for(i in seq_along(lCPT)){
    mm<-match(names(lout),colnames(lCPT[[i]]))
    mm<-colnames(lCPT[[i]])[mm[!is.na(mm)]]
    if(length(mm)>0){
     for(im in seq_along(mm)){
       lCPT[[i]]<- lCPT[[i]][!lCPT[[i]][, mm[im]] %in% lout[[mm[im]]],]
     }
    }

   }

   bn$CPTs <- lCPT
   return(bn)
 }
