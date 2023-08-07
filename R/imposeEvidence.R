#' imposeEvidence: a fuction for imposing evidence in the bayesian network.
#'
#' @param bn A bayesian network (output of buildBN function).
#' @import paramlink
#' @import igraph
#' @export
#' @return A preprocessed bayesian network.

imposeEvidence<-function(bn){
   G_SIZE_E     <- 8
   G_SIZE_X     <- 5
   G_SHAPE_E    <- "square"
   G_LABEL_DIST <- 0.4
   G_LABEL_CEX  <- 0.6
   G_ARROW_SIZE <- 0.6
   G_CURVED     <- 0.2

   E<-bn$E
   g<-bn$DAG
   EE <- list()
   isfounder<- which(unlist(lapply(strsplit(names(E),"_"),function(x){x[1]}))%in%bn$ped$founders)

   for(i in seq_along(E)){
     evtables <- paste0(names(E)[i],c("_p","_m"))
     homo     <- E[[i]][1]==E[[i]][2]

     if(!homo & !i%in%isfounder){
      newTableName<-paste0(names(E)[i],c("_pm"))
      cpt<-data.frame(rbind(c(E[[i]][1],E[[i]][2],paste(E[[i]][1:2],collapse="_"),1),
                            c(E[[i]][2],E[[i]][1],paste(E[[i]][2:1],collapse="_"),1)),stringsAsFactors=FALSE)
      colnames(cpt)<-c( paste0(names(E)[i],c("_p","_m")), newTableName,"prob")
      cpt[,"prob"]<-as.numeric(cpt[,"prob"])

      bn$CPTs[[newTableName]]<-cpt

      g<-set_vertex_attr(g,"variableType",evtables[1],"X")
      g<-set_vertex_attr(g,"variableType",evtables[2],"X")
      g<-set_vertex_attr(g,"size",evtables[1],G_SIZE_X)
      g<-set_vertex_attr(g,"size",evtables[2],G_SIZE_X)
      ccolor<-get.vertex.attribute(g,"color",evtables[1])
      g<-add.vertices(g,1,name=newTableName,locus=names(E)[i],pm="pm",nodeType="Gpm",
                      variableType="E",shape=G_SHAPE_E,size=G_SIZE_E,frame.color="white",
                      label.dist=G_LABEL_DIST,label.cex=G_LABEL_CEX,color=ccolor)
      g<-add.edges(g,c(evtables[1],newTableName, evtables[2],newTableName),
                     curved=G_CURVED,arrow.size=G_ARROW_SIZE)

      EE[[newTableName]]<-as.character(cpt[,newTableName])

      evidencia <- evtables[1]
      nei  <- neighbors(g,evidencia,mode="out")
      inei <- which(nei$nodeType=="G")
      if(length(inei)>0){
       nei <- nei[inei]
       for(j in seq_along(nei)){
        target <- names(nei)[j]
        cpt    <- bn$CPTs[[target]]
        cpt <- cpt[ cpt[,evidencia]%in%E[[i]],]
        bn$CPTs[[target]]<-cpt
       }
      }

      for(iet in 1:2){
       cpt <- bn$CPTs[[ evtables[iet] ]]
       cpt <- cpt[ cpt[,evtables[iet]]%in%E[[i]],]
       bn$CPTs[[evtables[iet]]] <- cpt
      }

     }else{
      for(ipm in 1:2){
       EE[[evtables[ipm]]] <- E[[i]][ipm]

       evidencia <- evtables[ipm]
       nei  <- neighbors(g,evidencia,mode="out")
       inei <- which(nei$nodeType=="G")
       if(length(inei)>0){
        nei <- nei[inei]
        for(j in seq_along(nei)){
         target <- names(nei)[j]
         cpt    <- bn$CPTs[[target]]
         cpt <- cpt[ cpt[,evidencia]%in%E[[i]][ipm],]
         bn$CPTs[[target]]<-cpt
        }
       }

       cpt <- bn$CPTs[[ evtables[ipm] ]]
       cpt <- cpt[ cpt[,evtables[ipm]]%in%E[[i]][ipm],]
       bn$CPTs[[evtables[ipm]]] <- cpt

      }
     }
   }
   bn$DAG <- g
   bn$clusters  <- clusters(g)
   bn$E   <- EE
   return(bn)

}
