#' getGenotypeTables: a function for obtaining genotypetables after variable elimination and using available genetic evidence.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param resQ List of CPTs.
#' @param geno data.frame with genotypes.
#' @param lqp list of individuals genotypes.
#' @import paramlink
#' @import igraph
#' @export
#' @return A dataframe with genotype probabilities.

getGenotypeTables<-function(bn,resQ,geno=NULL,lqp=NULL){

 if(is.null(lqp) & !is.null(geno)){
  lqp <- list()
  for(iqp in 1:nrow(geno)){
   saux    <- paste(rownames(geno)[iqp],colnames(geno),sep="_")
   alleles <- unlist(strsplit(geno[iqp,],"/"))
   for(i in seq_along(saux)){
    lqp[[saux[i]]]<- unname(alleles[(2*(i-1)+1):(2*i)])
   }
  }
 }
 if(is.null(lqp)){
  warning("Either lqp or geno should be specifid.\n")
  return()
 }

 rresQ<-resQ

 count<-1
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
      iprob <- grep("prob",colnames(S))
      if(length(iprob)>1){
       prob<-apply(S[,iprob,drop=FALSE],1,prod)
       S <- S[,-iprob]
       S <- cbind(S,prob=prob)
      }
      rresQ[[saux]] <- S
 }

 rresQ<-lapply(rresQ,function(x){
                icolpm <- grep("_pm",colnames(x))
                return(x[,-icolpm])
             })

 return(rresQ)
}
