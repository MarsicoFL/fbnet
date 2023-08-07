#' reportPQ: a function for calculating the probability of specified genotypes in a pedigree.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param resQ List of CPTs.
#' @param geno data.frame with genotypes.
#' @import paramlink
#' @import igraph
#' @export
#' @return A dataframe with genotype probabilities.

reportPQ<-function(bn,resQ,geno=NULL){
bn1 <- NULL
  if(is.null(geno)){
   pprod <- 1
   for(i in seq_along(resQ)){
    pprod <- pprod * sum(resQ[[i]][,"prob"])
   }
   pprod <- pprod *factorHeteroFounders(resQ,bn)
   LR1 <- pprod

   pprod <- 1
   for(i in seq_along(resQ)){
    pprod <- pprod * sum(resQ[[i]][,"prob"])
   }
   pprod <- pprod *factorHeteroFounders(resQ,bn)

   return(pprod)
  }

  lalleles <-apply(geno,1,function(x){strsplit(x,"/",fixed=TRUE)})
  saux<-c()
  for(i in seq_along(lalleles)){
   saux<-c(saux,paste(names(lalleles)[i],names(lalleles[[i]]),sep="_"))
  }
  lalleles<-unlist(lalleles,recursive=FALSE,use.names=FALSE)
  names(lalleles)<-saux


  haploNames <- gsub("_p","",gsub("_m","",names(resQ)))
  uhn <- unique(haploNames)

  resp<-c()
  for(i in seq_along(uhn)){

   ih   <- which(haploNames %in% uhn[i])
   ncol <- unlist(lapply(resQ[ih],ncol))

   personLoci <- strsplit(uhn[i],"__",fixed=TRUE)[[1]]
   nloci      <- length(personLoci)


    hap <- apply(matrix(unlist(lalleles[personLoci]),ncol=nloci),1,function(x){paste(x,collapse="_")})

    icols <- ncol[1] - (nloci:1)
    h1 <- apply(resQ[[ih[1]]][,icols,drop=FALSE],1,paste,collapse="_")
    icols <- ncol[2] - (nloci:1)
    h2 <- apply(resQ[[ih[2]]][,icols,drop=FALSE],1,paste,collapse="_")
    i1 <- c(which(h1%in%hap[1]),which(h2%in%hap[2]))
    i2 <- c(which(h1%in%hap[2]),which(h2%in%hap[1]))
    p<-0
    if(length(i1)==2) p <- p + resQ[[ih[1]]][i1[1],"prob"] * resQ[[ih[2]]][i1[2],"prob"]
    if(length(i2)==2) p <- p + resQ[[ih[1]]][i2[1],"prob"] * resQ[[ih[2]]][i2[2],"prob"]

#   }

    resp <- c(resp,p)
  }


  resp <- resp * factorHeteroFounders(resQ,bn1)

  res<-matrix(resp,ncol=1)
  rownames(res)<- uhn
  colnames(res)<-"likelihood"

  return(res)
}

