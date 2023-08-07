#' getLocusCPT: a function for obtaining the coditional probability table from a specific locus.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param locus Specified locus.
#' @param lumpingParameter Used for stepwise mutational model.
#' @param renorm If "row-wise" is selected, zero probability is assigned for transitions out of range. 
#' @import paramlink
#' @import igraph
#' @examples
#' pbn  <- initBN(toyped)
#' bnet <- buildBN(pbn,QP=3)
#' bn1  <- buildCPTs(bnet)
#' locCPT <- getLocusCPT(bn1,"M1")
#' @export
#' @return A bayesian network based on pedigree evidence and QP definition.

getLocusCPT<-function(bn,locus,lumpingParameter=NULL,renorm="row-wise"){ #,R=0.005,R2=5e-7,r=0.5){
  ilocus<-locus
  if(is.character(locus)){
    ilocus<-which(names(bn$alelFreq)==locus)
  }
  n<-length(bn$alelFreq[[ilocus]])
  a<-expand.grid(list(Ap=names(bn$alelFreq[[ilocus]]),Am=names(bn$alelFreq[[ilocus]]),S=1:2,node=names(bn$alelFreq[[ilocus]])),stringsAsFactors=FALSE)

  ih <- which(a[,"S"]==1)
  im <- which(a[,"S"]==2)

  b      <- rep(0,nrow(a))
  gender <- c("hombre","mujer")
  for(ig in seq_along(gender)){
    mmod <- names(bn$mmodel[[gender[ig]]])
    iag  <- which(a[,"S"]==ig)
    if(mmod=="none"){
      b[iag]<-apply(a[iag,],1,function(x){
        if(x["node"]==x[ig]) return(1)
        return(0) })
    }else{
      switch(mmod,
             equal={
               R <- bn$mmodel[[gender[ig]]][[1]]["R"]
               b[iag]<-apply(a[iag,],1,function(x){
                 if(x["node"]==x[ig]) return(1-R)
                 return(R/(n-1))
               })
             },
             stepwise={
               R  <- bn$mmodel[[gender[ig]]][[1]]["R"]
               R2 <- bn$mmodel[[gender[ig]]][[1]]["R2"]
               r  <- bn$mmodel[[gender[ig]]][[1]]["r"]

               nAlleles <- length(bn$alelFreq[[ilocus]])

               mij <- matrix(0,ncol=nAlleles,nrow=nAlleles)
               colnames(mij)<-rownames(mij)<-names(bn$alelFreq[[ilocus]])
               microVarGroup<-unlist(lapply(strsplit(colnames(mij),".",fixed=TRUE),function(x){if(length(x)==1) return(0);return(x[2])}))
               tt <- table(microVarGroup)
               s  <- sum(tt)-tt
               for(i in seq_along(mij[,1])){
                 mij[i,i] <- 1 - R - R2
                 imv   <- which(microVarGroup==microVarGroup[i])
                 step  <- as.numeric( colnames(mij)[imv]) - as.numeric( colnames(mij)[i])
                 knorm <- R/sum(r^abs(step[step!=0]))
                 for(j in 1:ncol(mij)){
                   if(j==i) next
                   if(microVarGroup[i]!=microVarGroup[j]){
                     mij[i,j] <- R2/s[microVarGroup[i]]
                   }else{
                     step <- as.numeric(colnames(mij)[j])- as.numeric( colnames(mij)[i])
                     mij[i,j] <-knorm*r^abs(step)
                   }
                 }
               }

               if(names(bn$mmode[[1]])!="none" & !is.null(lumpingParameter)){
                 minp <- 0
                 if(is.logical(lumpingParameter)){
                   rn <- cn <- as.numeric(colnames(mij))
                   pmin<-c()
                   for(i in seq_along(rn)){
                     j <- which(abs(cn-rn[i])==1)
                     pmin <- c(pmin,mij[i,j])
                   }
                   minp<-min(pmin)
                 }else{
                   if(lumpingParameter<1){
                     minp <- lumpingParameter
                   }
                 }
                 if(minp>0){  

                   mij <- apply(mij,1,function(x){
                     x[x<minp]<-0
                     return(x/sum(x))
                   })

                 }else{
                   rn <- cn <- as.numeric(colnames(mij))
                   maux<-mij

                   for(i in seq_along(rn)){
                     dif  <- rn[i]-cn
                     iout <- which(!(dif==round(dif) & abs(dif)<=lumpingParameter))
                     iin  <- which((dif==round(dif) & abs(dif)<=lumpingParameter))

                     maux[i,-iin]<-0

                     if(renorm=="row-wise"){
                       maux[i,]<-maux[i,]/sum(maux[i,])
                     }else{
                       auxNorm <- 1-maux[i,i]
                       maux[i,iin[!iin%in%i]]<-maux[i,iin[!iin%in%i]]/sum(maux[i,iin[!iin%in%i]])*auxNorm
                     }
                   }
                   mij <- maux

                 }
               }


               b[iag]<-apply(a[iag,],1,function(x){
                 saux<-as.character(x)[c(ig,length(x))]
                 return(mij[saux[1],saux[2]])
               })
             }
      )
    }
  }
  return(cbind(a,prob=b))
}

