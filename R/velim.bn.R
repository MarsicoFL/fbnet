#' velim.bn: a function for variable elimination in a bayesian network.
#'
#' @param bn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param ordMethod Selected ordering method between id, min_degree, min_fill and fixed.
#' @param orderElim Elimination order.
#' @param verbose Computation output.
#' @import paramlink
#' @import igraph
#' @examples
#' pbn  <- initBN(toyped)
#' bnet <- buildBN(pbn,QP=3)
#' bn1  <- buildCPTs(bnet)
#' resQ <- velim.bn(bn1,ordMethod="min_fill",verbose=FALSE)
#' @export
#' @return Variable elimination result.

velim.bn <- function(bn,ordMethod=c("id","min_degree","min_fill","fixed")[2],orderElim=NULL,verbose=FALSE){
  Q  <- bn$Q
  E  <- bn$E
  df <- NULL 
intersection <- function(x, y, ...){
     if (missing(...)) intersect(x, y)
 else intersect(x, intersection(y, ...))
 }


setOrdering<-function(bn,ordMethod,vars=NULL,orderElim=NULL){

 bn$ordMethod <- ordMethod 
 
 if(ordMethod=="id"){
  if(!is.null(vars)){
    bnordering<-names(bn$vars)[!names(bn$vars)%in%bn$Q] 
    bnordering<-bn[["ordering"]][bn[["ordering"]]%in%vars]
  }else{
   bnordering<-names(bn$vars)[!names(bn$vars)%in%bn$Q] 
  } 
 }
 
 if(ordMethod%in%c("min_degree","min_fill")){
   bnordering<-.minOrdering(bn,vars,ordMethod)
 }
 
 if(ordMethod=="fixed"){
  if(is.null(orderElim)) warning("ordMethod=fixed, but no orderElim provided.\n")
  if(!all(orderElim%in%names(bn$CPTs))) warning("check orderElim bariable names.\n")
  bnordering<-orderElim
 } 
 
 return(bnordering)  
}

.prodFactor.orig<-function(laux,verbose=FALSE){

 if(length(laux)==1) return(laux[[1]])


 #identifico las tablas que en realidad son un escalar
 iscalar <- which(unlist(lapply(laux,is.numeric)))
 kscalar <- 1
 if(length(iscalar)>0){
  kscalar <- prod(unlist(laux[iscalar]))
  laux<-laux[-iscalar]  # sigo procesando las tabals que no son escalares 
 } 
 
 
 
 if(length(laux)==0) return(kscalar)
 
 

 
 z <- lapply(laux,function(x){return(colnames(x))})
 z <- unique(c(names(z),unlist(z)))
 z <- z[z!="prob"]
 lzvars <- list()
 for(zz in z){
  varstab <- colnames(laux[[zz]])
  varstab <- varstab[!varstab%in%"prob"]
  for(iv in seq_along(varstab)){
   lzvars[[varstab[iv]]] <- sort(unique(laux[[zz]][,varstab[iv]]))
  } 
 }
 z<-expand.grid(lzvars,stringsAsFactors=FALSE)
 f<-rep(1,nrow(z))
 
 for(iz in seq_along(z[,1])){
  if(verbose)cat(iz,"/",nrow(z),"\n")
  #por cada instanciacion recorro tablas
  for(i in 1:length(laux)){
    
    #identifico que variables instanciadas aparecen en la tabla
    zcols <- colnames(z)[colnames(z)%in%colnames(laux[[i]])]
    
    #identifico la fila de la tabla
    irow<-which(apply(laux[[i]],1,function(x){ all(x[zcols]==z[iz,zcols])}))
    if(length(irow)>0){
      f[iz] <- f[iz] * laux[[i]][irow,"prob"]
    }else{
      f[iz] <- 0    ## ESTA BIEN ESTO, NO?
      break
    }
    
  }
  
 }
 df <- data.frame(z,prob=f*kscalar)
 colnames(df)<-c(colnames(z),"prob")
 return(df)
}


T1 <- expand.grid(list(A=1:4,B=1:2,C=1:3))
T2 <- expand.grid(list(B=1:5,C=1:2,D=1:3))
T3 <- expand.grid(list(A=c(2,3,5),B=1:3,C=1:3,E=1:2))
T1<-cbind(T1,prob=stats::runif(nrow(T1)))
T2<-cbind(T2,prob=stats::runif(nrow(T2)))
T3<-cbind(T3,prob=stats::runif(nrow(T3)))
laux<-list(T1=T1,T2=T2,T3=T3)
 
.prodFactor.brute.force.incomplete<-function(laux){

 if(length(laux)==1) return(laux[[1]])


 iscalar <- which(unlist(lapply(laux,function(x){ 
                                 res<-nrow(x)
                                 if(is.null(res)) return(TRUE)
                                 if(res==1) return(TRUE)
                                 return(FALSE)
                               })))
 kscalar <- 1
 if(length(iscalar)>0){
  kscalar <- prod(unlist(laux[iscalar]))
  laux<-laux[-iscalar]  # sigo procesando las tabals que no son escalares 
 } 
 
 
 
 if(length(laux)==0) return(kscalar)
 
 z <- lapply(laux,function(x){return(colnames(x)[!colnames(x)%in%"prob"])})
 
 { 
  commonVarNames<-unlist(lapply(reverseSplit(z),function(x){length(x)==length(laux)}))
  commonVarNames<-names(commonVarNames)[commonVarNames & names(commonVarNames)!="prob"]
  lCommonCombo<-lapply(laux,function(x){ 
                 return(apply(x[,commonVarNames],1,paste0,collapse="."))
               })
 
  ll<-reverseSplit(lCommonCombo)
  commonVars<-names(which(unlist(lapply(ll,function(x){length(unique(x))==length(laux)}))))
  if(length(commonVars)==0) return(0) #no hay combinacion comun minima
  
  lKeep<-lapply(lCommonCombo,function(x){
           return(which(x%in%commonVars))
          })
          
  for(i in seq_along(laux)){
   laux[[i]]<-laux[[i]][lKeep[[i]],]
  }
 } 
 
 
 { 
  nonCommon <- unique(unlist(lapply(laux,function(x){return(colnames(x)[!colnames(x)%in%commonVarNames])})))
  nonCommon <- nonCommon[!nonCommon%in%"prob"]
  lvars     <- lapply(seq_along(nonCommon),function(x){return(-1)})
  names(lvars)<-nonCommon
  for(l in seq_along(laux)){
    x <- laux[[l]]
    ccol <- colnames(x)[!colnames(x)%in%commonVarNames]
    ccol <- ccol[!ccol%in%"prob"]
    if(length(ccol)==0) next
    for(i in seq_along(ccol)){
      if(min(lvars[[ccol[i]]])==-1){
       lvars[[ccol[i]]]<-unique(x[,ccol[i]])
      }else{
       lvars[[ccol[i]]]<-unique(intersect(lvars[[ccol[i]]],x[,ccol[i]]))
      } 
    }
 }
 
  lres<-lapply(laux,function(x){
            ccol <- colnames(x)[colnames(x)%in%names(lvars)]
            ii <- rep(TRUE,nrow(x))
            for(i in seq_along(ccol)){
             ii<- ii  & x[,ccol[i]]%in%lvars[[ccol[i]]]
            }
            return(x[ii,])
          })
  }        
  
 
 
 
 pprod <- lres[[1]]
 if(length(lres)>1){
  for(i in 2:length(lres)){
   ccol <- intersect(colnames(pprod),colnames(lres[[i]]))
   ccol <- ccol[!ccol%in%"prob"]
   pprod <- merge(pprod,lres[[i]],by=ccol)
  }
 } 
 
 t(apply(pprod,1,function(x){
            iprob <- grep("prob",names(x))
            return( c(x[-iprob],prob=prod(x[iprob])) )
          }))
  
  
 return(df)
}

.prodFactor<-function(laux){

 if(length(laux)==1) return(laux[[1]])
 pprod <- laux[[1]]
 
 if(length(laux)>1){
  for(i in 2:length(laux)){
   ccol <- intersect(colnames(pprod),colnames(laux[[i]]))
   ccol <- ccol[!ccol%in%"prob"]
   pprod <- suppressWarnings(merge(pprod,laux[[i]],by=ccol))
  }
  iprob <- grep("prob",colnames(pprod))
  df <-cbind(pprod[,-iprob],prob=apply(pprod[,iprob],1,prod)) 
 }else{
  df <- pprod
 }
 
 return(df)
}

.sumFactor<-function(cpt,Z){
 X <- colnames(cpt)
 Y <- X[!X%in%c(Z,"prob")]
 if(length(Y)==0) return(1)  #Z=Y!!
 

 if(length(Y)>1){
  f<-stats::aggregate(cpt[,"prob"],by=as.list(cpt[,Y]),sum)
 }else{ 
  f<-stats::aggregate(cpt[,"prob"],by=list(cpt[,Y]),sum)
 } 

 colnames(f)<-c(Y,"prob")
 return(f)
}

.minDegreeOrdering<-function(bn,vars=NULL){
 g <- make_empty_graph(directed=FALSE)
 
 if(is.null(vars)){ 
  lcpt <- bn$CPTs
 }else{
  lcpt <- bn$CPTs[vars]
 }
 
 for(i in seq_along(lcpt)){
      vcond<-colnames(lcpt[[i]])
      if(!is.null(vcond)){
       vcond<-vcond[!vcond%in%c("prob",bn$Q)]
       if(length(vcond)>0){
         gg <- make_full_graph(length(vcond))
         V(gg)$name <- vcond
         g <- g + gg
       }
      }
  }    
  allvars <- V(g)$name
  
  ord <- c()
  for(i in 1:(vcount(g)-1)){
    imin<-which.min(degree(g))
    ord <- c(ord,names(imin))
    nn  <- neighbors(g,imin)
    if(length(nn)>1){  #enlazo vecinos del que voy a eliminar
      gg <- make_full_graph(length(nn))
      V(gg)$name <- names(nn)
      g <- g + gg
    }
    g <- delete_vertices(g,names(imin))
    #plot(g)
  }
  ord <- c(ord,setdiff(allvars,ord))
  
  return(ord)
}

.minOrdering<-function(bn,vars=NULL,method=c("min_degree","min_fill")[1]){
 g <- make_empty_graph(directed=FALSE)
 
 if(is.null(vars)){ 
  lcpt <- bn$CPTs
 }else{
  lcpt <- bn$CPTs[vars]
 }
 
 for(i in seq_along(lcpt)){
      vcond<-colnames(lcpt[[i]])
      if(!is.null(vcond)){
       vcond<-vcond[!vcond%in%c("prob",bn$Q)]
       if(length(vcond)>0){
         gg <- make_full_graph(length(vcond))
         V(gg)$name <- vcond
         g <- g + gg
       }
      }
  }    
  allvars <- V(g)$name
  
  ord <- c()
  for(i in 1:(vcount(g)-1)){
    if(method=="fill"){     
     tc <- transitivity(g,"local") 
     tc[is.na(tc)] <- 1000         
     dg <- degree(g)
     numVecNoConectados<-(1-tc)*dg*(dg-1)/2
     imin <- which.min(numVecNoConectados)
    }else{
     imin<-which.min(degree(g))
    } 
    ord <- c(ord,names(imin))
    nn  <- neighbors(g,imin)
    if(length(nn)>1){  #enlazo vecinos del que voy a eliminar
      gg <- make_full_graph(length(nn))
      V(gg)$name <- names(nn)
      g <- g + gg
    }
    g <- delete_vertices(g,names(imin))
  }
  ord <- c(ord,setdiff(allvars,ord))
  
  return(ord)
}


 if(length(Q)==0){warning("por ahora falta chequear implementacion a nivel sistema. igual no hace falta para genis\n")}
 
 lf<-list()
 for(ic in seq_along(1:bn$clusters$no)){
  vars <- names(bn$clusters$membership)[bn$clusters$membership==ic]
  vars <- vars[vars%in%names(bn$CPTs)] 
  

  S  <- bn[["CPTs"]][vars]
 
  S <- lapply(S,function(x){
              return(x[x$prob>0,,drop=FALSE])
             })
             
  oord <- setOrdering(bn,ordMethod,vars,orderElim)
  oord <- oord[oord%in%vars] #hace falta cuando se pasa un ordering desde afuera
  
  for(iord in seq_along(oord)){
   ik   <- which(unlist(lapply(names(S),function(x){return(oord[iord]%in%c(x,colnames(S[[x]])))})))
   faux <- .prodFactor(S[ik])

   if(verbose){
     cat(iord,oord[iord],paste("factor:",paste(names(S)[ik],collapse=":")))
   }
   
   if(nrow(faux)>1){
    fi   <- .sumFactor(faux,oord[iord])
    S[ik] <- NULL
    S[[paste0("f",iord)]]<-fi
    
    if(verbose) cat(paste0(" -> f",iord," width:",nrow(fi)))
   }
   if(verbose) cat("\n")
  }
  
  lf[[ic]] <- .prodFactor(S)
 }
   
 if(length(Q)==0){
   if(FALSE){
    res <- prod(unlist(lapply(lf,function(x){x$prob})))
   
    fnds<-paste0(bn$ped$founders,"_")
    num_hetero<-0
    for(i in seq_along(fnds)){
     ii<-grep(fnds[i],names(bn$E))
     if(length(ii)%%2!=0) warning(paste("Ojo...nro impar de alelos para el founder:",fnds[i])) 
     num_hetero<-num_hetero+sum(apply(matrix(bn$E[ii],byrow=TRUE,ncol=2),1,function(x){x[1]!=x[2]}))
    }
    res <-res*2^num_hetero
   }
  
  sysName  <- p <- num_hetero <- c()
  fnds<-paste0(bn$ped$founders,"_") 
  for(i in seq_along(lf)){
    este<-strsplit(names(lf[[i]])[1],"_",fixed=TRUE)[[1]][2]

     i1<-grep(este,names(bn$E))
     i2<-c()
     for(j in seq_along(fnds)) i2<-c(i2,grep(fnds[j],names(bn$E)))
     ii<-intersect(i1,i2)
     m<-as.data.frame(matrix(bn$E[ii],byrow=TRUE,ncol=2),drop=FALSE)
     
     
     num_hetero<-c(num_hetero,sum(apply(m,1,function(x){x[1]!=x[2]})))
    sysName<-c(sysName,este)
    p <- c(p,lf[[i]][,"prob"])
  }

  res<-cbind(by(cbind(p,num_hetero),sysName,function(x){return(prod(x[,1])*2^(sum(x[,2])/2))}))
  colnames(res)<-c("likelihood")
  

  
 }else{

  names(lf) <- unlist(lapply(lf,function(x){paste(colnames(x)[colnames(x)%in%Q],collapse="__")}))
  
  res<-lf
 }
 res<-lapply(res,function(x){
              x<-x[,colnames(x)%in%c(bn$Q,"prob")]
            })
 return(res)
}
