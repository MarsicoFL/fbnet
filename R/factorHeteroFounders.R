#' factorHeteroFounders: a function for multiplying probabilities in case of heterocigote founders.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param rresQ List of CPTs.
#' @import paramlink
#' @import igraph
#' @export
#' @return A dataframe with genotype probabilities.

factorHeteroFounders<-function(rresQ,bn){

  fnds<-paste0(bn$ped$founders,"_")

  haploNames <- sub("^\\d+_","",names(rresQ),fixed=FALSE)
  haploNames <- sub("__\\d+_","__",haploNames,fixed=FALSE)
  hn <- unique(gsub("_p","",gsub("_m","",haploNames)))

  nHetero<-rep(0,length(hn))
  names(nHetero)<-hn
  for(i in seq_along(fnds)){
   haploNames <- sub("^\\d+_",fnds[i],names(rresQ),fixed=FALSE)
   haploNames <- sub("__\\d+_",paste0("__",fnds[i]),haploNames,fixed=FALSE)

   lhaplo <- strsplit(unique(gsub("_p","",gsub("_m","",haploNames))),"__")
   for(j in seq_along(lhaplo)){
    aacum<-c()
    for(k in seq_along(lhaplo[[j]])){
      ia       <- grep(lhaplo[[j]][k],names(bn$E))
      aacum    <- paste(aacum,bn$E[ia])
    }
    isHetero <- length(unique(aacum))>1
    nHetero[j] <- nHetero[j] + as.numeric(isHetero)
   }
  }
  res <- 2^nHetero
  return(res)
}

