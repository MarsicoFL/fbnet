#' reportLR: a function for calculating the LRs of specified genotypes in a pedigree.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param resQ List of CPTs.
#' @param geno data.frame with genotypes.
#' @import paramlink
#' @import igraph
#' @export
#' @return A dataframe with LRs.

reportLR<-function(bn,resQ,geno=NULL){


 rresQ<-resQ

 lqp <- list()
 for(iqp in 1:nrow(geno)){
  saux    <- paste(rownames(geno)[iqp],colnames(geno),sep="_")
  alleles <- unlist(strsplit(geno[iqp,],"/"))
  for(i in seq_along(saux)){
   lqp[[saux[i]]]<- unname(alleles[(2*(i-1)+1):(2*i)])
  }
 }

 rresQ <- getGenotypeTables(bn,resQ,lqp=lqp)


if(FALSE){
 count<-rmp<-1
 for(i in seq_along(lqp)){
      newTableName<-paste0(names(lqp)[i],c("_pm"))
      if(lqp[[i]][1]!=lqp[[i]][2]){
       cpt<-data.frame(rbind(c(lqp[[i]][1],lqp[[i]][2],paste(lqp[[i]][1:2],collapse="_")),
                             c(lqp[[i]][2],lqp[[i]][1],paste(lqp[[i]][2:1],collapse="_"))),stringsAsFactors=FALSE)
      }else{
       cpt<-data.frame(rbind(c(lqp[[i]][1],lqp[[i]][2],paste(lqp[[i]][1:2],collapse="_"))),stringsAsFactors=FALSE)
      }
      colnames(cpt)<-c( paste0(names(lqp)[i],c("_p","_m")), newTableName)

      pm2table<-unlist(lapply(rresQ,function(x){return(any(colnames(x)%in%colnames(cpt)))}))
      pm2table<-names(pm2table)[pm2table]

      S <- cpt
      saux <-c()
      for(j in seq_along(pm2table)){
          bycol <- intersect(colnames(rresQ[[pm2table[j]]]),colnames(S))
          bycol <- bycol[!bycol%in%"prob"]
          S <- merge(S,rresQ[[pm2table[j]]],by=bycol)
          rresQ[[pm2table[j]]]<-NULL
          count<-count+1
          if(j==1){
           saux =pm2table[j]
          }else{
           saux <- paste(saux,pm2table[j],sep="__")
          }
      }
      rresQ[[saux]] <- S

      sys <- unlist(lapply(strsplit(names(lqp)[i],"_"),function(x){return(x[2])}))

      rmp <- rmp * prod(bn$alelFreq[[sys]][lqp[[i]]]) * 2^(lqp[[i]][1]!=lqp[[i]][2])
  }

}

 fhetero<-factorHeteroFounders(resQ,bn)

 systemTable1<-lapply(rresQ,function(x){
                       xx<-strsplit(colnames(x),"_",fixed=TRUE)
                       xx<-xx[lapply(xx,length)==3]
                       saux<-sort(unique(unlist(lapply(xx,function(y){return(y[2])}))))
                       return(saux)})

 names(systemTable1)<-unlist(lapply(systemTable1,function(x){
                                saux<-paste(sort(x),collapse="_")
                              }))

 lL1 <- lapply(rresQ,function(x){sum(x[,"prob"])})



 lL0<-list()
 for(i in seq_along(systemTable1)){
   rmp<-1
   for(j in 1:length(systemTable1[[i]])){
    isysq <- grep(systemTable1[[i]][j],names(lqp))
    for(iq in seq_along(isysq)){
     isq <- isysq[iq]
     rmp <- rmp * prod(bn$alelFreq[[systemTable1[[i]][j]]][lqp[[isq]]]) * 2^(lqp[[isq]][1]!=lqp[[isq]][2])
    }
   }

   isys0 <- which(unlist(lapply(strsplit(names(resQ),"_"),function(x){x[2]}))%in%systemTable1[[i]])
    lL0[[names(systemTable1)[i]]] <- prod(unlist(lapply(resQ[isys0],function(x){
                                                   return(sum(x[,"prob"]))})))*rmp
 }

 for(i in seq_along(lL1)){
  lL1[[i]] <- lL1[[i]]*fhetero[i]
  lL0[[i]] <- lL0[[i]]*fhetero[i]
 }


 LR <- prod(unlist(lL1))/prod(unlist(lL0))

 df <- data.frame(L0=unlist(lL0),L1=unlist(lL1),LR=unlist(lL1)/unlist(lL0))
 rownames(df)<-names(systemTable1)

 return(list(systemLR=df,LR=LR))

}
