#' buildCPTs: a function for building conditional probability tables based on pedigree bayesian network.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param bNodePrunning Standard pruning.
#' @param bStateRemoval State based pruning.
#' @param bStateRemoval2 State based pruning (model 2).
#' @param lumpingParameter Used for stepwise mutational model.
#' @param renorm If "row-wise" is selected, zero probability is assigned for transitions out of range. 
#' @param verbose Computations output. 
#' @import paramlink
#' @import igraph
#' @examples
#' pbn  <- initBN(toyped)
#' bnet <- buildBN(pbn,QP=3)
#' bn1  <- buildCPTs(bnet)
#' @export
#' @return A bayesian network based on pedigree evidence and QP definition.

buildCPTs<-function(bn,bNodePrunning=TRUE,bStateRemoval=TRUE,bStateRemoval2=TRUE,lumpingParameter=NULL,renorm="row-wise",verbose=FALSE){
  debug<-verbose

  lAllelesCPTs <-list()
  for(i in seq_along(bn$alelFreq)){
   locus <- names(bn$alelFreq)[i]
   paux <- getLocusCPT(bn,locus,lumpingParameter,renorm)

   lAllelesCPTs[[locus]] <- paux
  }



  lSelectorCPTs<-list()
  for(i in seq_along(bn$markerLinkage[,1])){
    if(bn$markerLinkage[i,"recomb"]!=0.5){
     r <- bn$markerLinkage[i,"recomb"]
     m <- matrix(c(r,1-r,1-r,r),2,2)
     colnames(m)<-paste("AS",c(1:2),sep="_") 
     rownames(m)<-paste("S",c(1:2),sep="_")
     cpt<-expand.grid(list(AS=c(1,2),S=c(1,2)),stringsAsFactors=FALSE)
     cpt<-data.frame(cpt,prob=apply(cpt,1,function(x){m[x[1],x[2]]}),stringsAsFactors=FALSE)

     saux<-paste(rownames(bn$markerLinkage)[i],bn$markerLinkage[i,"linkedTo"],sep="_")
     lSelectorCPTs[[saux]]<-cpt
    }
   }

  d     <- degree(bn$DAG,mode="in")
  lCPTs <- list()
  for(i in seq_along(d)){
   variableName<-names(d)[i]
   if(d[i]==0){                                    
    if(V(bn$DAG)[variableName]$nodeType=="G"){     
     a<-bn$alelFreq[[V(bn$DAG)$locus[i]]]
     cpt          <- data.frame(names(a),a,stringsAsFactors=FALSE)
     colnames(cpt)<- c(variableName,"prob")
     lCPTs[[variableName]]<-cpt
    }else{                                         
     cpt<-data.frame(c(1,2),c(.5,.5),stringsAsFactors=FALSE)
     colnames(cpt)<-c(variableName,"prob")
     lCPTs[[variableName]]<-cpt
    }
   }else{
    nei<-neighbors(bn$DAG,V(bn$DAG)[variableName],mode="in")
    if(V(bn$DAG)[variableName]$nodeType=="G"){
     ap <- nei$name[which(nei$nodeType=="G" & nei$pm=="p")] 
     am <- nei$name[which(nei$nodeType=="G" & nei$pm=="m")] 
     s  <- nei$name[which(nei$nodeType=="S")]               

     cpt <- lAllelesCPTs[[V(bn$DAG)[variableName]$locus]]
     colnames(cpt)[1:4]<-c(ap,am,s,variableName)

     lCPTs[[variableName]]<-cpt
    }else{                                         
     aux<-paste(V(bn$DAG)[variableName]$locus,nei$locus,sep="_")
     cpt<-lSelectorCPTs[[aux]]
     colnames(cpt)<-c(nei$name,variableName,"prob")
     lCPTs[[variableName]]<-cpt
    }
   }
  }
  bn$CPTs <- lCPTs
  if(debug){
   foo<-function(x){
      a<-unlist(lapply(x$CPTs,function(xx){sum(xx$prob>0)}))
      res<-rep(0,length(lCPTs));names(res)<-names(lCPTs)
      res[names(a)]<-a
      return(res)
    }
    mfoo<-c()
    mfoo<-cbind(mfoo,orig=foo(bn))
  }




  if(bNodePrunning){
   bn<-pruneNodes(bn)
   if(debug) mfoo<-cbind(mfoo,nodePrunning=foo(bn))
  }

  bn<-imposeEvidence(bn)
  if(debug){
    mfoo1<-mfoo
    mfoo2<-foo(bn)
    if(length(mfoo2)>nrow(mfoo1)){
     nc<-ncol(mfoo1)
     nr<-length(mfoo2)-nrow(mfoo1)
     mfoo<-rbind(mfoo1,matrix(rep(0,nc*nr),ncol=nc,nrow=nr))
    }
    mfoo<-cbind(mfoo,imposeEvidence=foo(bn))
    rownames(mfoo)<-names(mfoo2)
  }



  if(bStateRemoval){
    bn <- stateRemoval(bn)
    if(debug) mfoo<-cbind(mfoo,stateRemoval=foo(bn))
  }
  if(bStateRemoval2){
   bn <- stateRemoval2(bn)
   if(debug) mfoo<-cbind(mfoo,stateRemoval2=foo(bn))
  }



  lvars <- list()
  for(l in seq_along(bn$CPTs)){
    x <- bn$CPTs[[l]]
    for(i in (1:(ncol(x)-1))){
      lvars[[colnames(x)[i]]]<-unique(c(lvars[[colnames(x)[i]]],x[,colnames(x)[i]]))
    }
  }
  bn$vars <- lvars
  if(debug)message(mfoo)
  return(bn)
 }
