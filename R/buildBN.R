#' buildBN: a function for building the bayesian network.
#'
#' @param pbn A bayesian network for pedigree object with information of the genotyped members. The ped object must be in Familias format.
#' @param QP Query Persons Ids
#' @import paramlink
#' @import igraph
#' @import grDevices
#' @examples
#' pbn  <- initBN(toyped)
#' bnet <- buildBN(pbn,QP=3)
#' @export
#' @return A bayesian network based on pedigree evidence and QP definition.

buildBN<-function(pbn,QP){
  G_SIZE_E     <- 8
  G_SIZE_X     <- 5
  G_SHAPE_E    <- "square"
  G_LABEL_DIST <- 0.4
  G_LABEL_CEX  <- 0.6
  G_ARROW_SIZE <- 0.6
  G_CURVED     <- 0.2

  ped <- pbn$ped
  markerEvidence=pbn$markerEvidence
  mutationModel=pbn$mmodel
  markerLinkage=pbn$markerLinkage
  alelFreq=pbn$alelFreq   
  
  
  bn<-list()
  
  res<-c()
  pped <- ped$ped
  for(iid in seq_along(pped[,1])){
   for(ilocus in seq_along(alelFreq)){
    locusName   <- names(alelFreq)[ilocus]
    linkedTo<-NULL
    if(!is.null(markerLinkage)){
     if(markerLinkage[locusName,"linkedTo"]!=""){
      linkedTo <- markerLinkage[locusName,"linkedTo"]
      recomb   <- markerLinkage[locusName,"recomb"]
     } 
    }
    
    for(pm in c("p","m")){
      nodeName <- paste(pped[iid,"ID"],locusName,pm,sep="_")
      
      if(pm=="p" & pped[iid,"FID"]!=0){
       for(apm in c("p","m")){
        nodeNameFrom <- paste(pped[iid,"FID"],locusName,apm,sep="_")
        res <- c(res,nodeNameFrom,nodeName)
       } 
       nodeNameSelector <- paste(pped[iid,"ID"],locusName,pm,"S",sep="_")
       res<-c(res,nodeNameSelector,nodeName)
       
      
       if(!is.null(linkedTo)){
        nodeNameFrom <- paste(pped[iid,"ID"],linkedTo,pm,"S",sep="_")
        res <- c(res, nodeNameFrom, nodeNameSelector)
       }       
      }
      
      if(pm=="m" & pped[iid,"MID"]!=0){
       for(apm in c("p","m")){
        nodeNameFrom <- paste(pped[iid,"MID"],locusName,apm,sep="_")
        res <- c(res,nodeNameFrom,nodeName)
       } 
       nodeNameSelector <- paste(pped[iid,"ID"],locusName,pm,"S",sep="_")
       res<-c(res,nodeNameSelector,nodeName)
      
       if(!is.null(linkedTo)){
        nodeNameFrom <- paste(pped[iid,"ID"],linkedTo,pm,"S",sep="_")
        res <- c(res, nodeNameFrom, nodeNameSelector)
       }        
      }
      
      

    }   
   }  

  }
 
  nodes <- sort(unique(res))
  df<-data.frame(stringsAsFactors=FALSE)
  for(i in seq_along(nodes)){
   aux<-strsplit(nodes[i],"_",fixed=TRUE)[[1]]
   if(length(aux)==3) aux<-c(aux,"G")
   aaux<-data.frame(nodeName=nodes[i],
                    id=aux[1],
                    locus=aux[2],
                    pm=aux[3],
                    nodeType=aux[4],stringsAsFactors=FALSE)                       
   df<-rbind(df,aaux)
  }
  rownames(df)<-df[,1]
 
  Q<-df[df[,"id"]%in%QP & df[,"nodeType"]=="G","nodeName"]
  

  
  g <- graph_from_data_frame(matrix(res,ncol=2,byrow=TRUE),vertices=df)
  
  
  E<-list()
  for(j in seq(1,ncol(markerEvidence),2)){
   for(i in seq_along(markerEvidence[,1])){
    e<-unname(markerEvidence[i,c(j,j+1)])
    if(all(is.na(e))) next
    saux <- paste(rownames(markerEvidence)[i],colnames(markerEvidence)[j],sep="_")
    saux <- gsub("_m","",gsub("_p","",saux))
    E[[saux]]<-e
   }
  }    
  
  bn$ped            <- ped
  bn$markerEvidence <- markerEvidence
  bn$mmodel <- mutationModel
  bn$markerLinkage  <- markerLinkage
  bn$alelFreq       <- alelFreq
  bn$E              <- E
  bn$Q              <- Q
  
  
  V(g)$variableType           <- "X"
  V(g)[paste(names(E),c("p","m"),sep="_")]$variableType <- "E"
  V(g)[Q]$variableType        <- "Q"
  
  
  V(g)$shape       <- ifelse(df[V(g)$name,"nodeType"]=="G",G_SHAPE_E,"circle")
  
  V(g)$size                         <- G_SIZE_E
  V(g)[V(g)$nodeType=="S"]$size     <- G_SIZE_X
  V(g)[V(g)$variableType=="X"]$size <- G_SIZE_X
  
  
  V(g)$frame.color    <- "white"
  V(g)[Q]$frame.color <- "black"
  
  V(g)$label.dist    <- G_LABEL_DIST
  V(g)$label.cex     <- G_LABEL_CEX
  
  u<-unique(V(g)$id)
  V(g)$color       <- grDevices::rainbow(length(u),alpha=0.5)[match(V(g)$id,u)]

  E(g)$curved      <- G_CURVED
  E(g)$arrow.size  <- G_ARROW_SIZE
  
  bn$DAG <- g

  return(bn)
 } 
