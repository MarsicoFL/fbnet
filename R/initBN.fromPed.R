#' initBN.fromPed: a function to initialize the bayesian network.
#'
#' @param bplotped An alternative ped object to be compared. 
#' @param ped A ped object in Familias format. 
#' @import paramlink
#' @import graphics
#' @export
#' @return A bayesian network.

initBN.fromPed<-function(ped,bplotped){
  lLociFreq <-  NULL

persons <- as.character(ped$orig.ids)
pid     <- ped$ped[,"FID"];pid[pid==0]<-NA
mid     <- ped$ped[,"MID"];mid[mid==0]<-NA
sex     <- ifelse(ped$ped[,"SEX"]==1,"male","female")
QP      <- ped$pedigree[ped$pedigree[,"AFF"]==2,"ID"]
systems <- unlist(lapply(ped$markerdata,attr,"name"))
linkageR <- rep(0.5,length(systems))

bSimuData  <- FALSE
knownIds <- ped$available

ped1 <- FamiliasPedigree(id=persons,dadid=pid,momid=mid,sex=sex)



myloci <- list()
for(i in seq_along(systems)){
 if(!systems[i]%in%names(lLociFreq)) next
 freqs <- lLociFreq[[systems[i]]][,"freq"]
 anames<- lLociFreq[[systems[i]]][,"Alelo"]
 anames[anames=="0"]<-"zero"
 if(sum(freqs)>1) freqs<-freqs/sum(freqs)
 if(1-sum(freqs)>1e-4){
  freqs<-c(freqs,1-freqs)
  anames<-c(anames,"ExtraAlelle")
 }
 locus<-FamiliasLocus(frequencies=freqs,
                      allelenames=anames,
                      name=systems[i])
 myloci[[systems[i]]]<-locus
}

markerLinkage<- data.frame(linkedTo=rep("",length(systems)),
                           recomb=rep(0.5,length(systems)),stringsAsFactors=FALSE)
rownames(markerLinkage)<-systems
for(i in seq_along(linkageR)){
 if(linkageR[i]!=0.5 & i>1){
  markerLinkage[i,"linkedTo"] <- rownames(markerLinkage)[i-1]
  markerLinkage[i,"recomb"]   <- linkageR[i]
 }
}


auxped <- paramlink::Familias2linkdat(ped1,NULL,myloci)
auxped$markerdata <- ped$markerdata
auxped$available  <- ped$available
auxped$nMark      <- ped$nMark

if(!is.null(QP)){


 auxped$pedigree[auxped$orig.ids%in%QP,"AFF"]<-2
}

if(bplotped){
	graphics::plot(auxped,mark=1:min(3,length(myloci)))
 graphics::mtext(paste("H1:",length(myloci),"markers"))
}


 alelFreq <- lapply(auxped$markerdata,function(x){
                         res<-attr(x,"afreq");
                         names(res)<-attr(x,"alleles")
                         return(res)
                         })
 names(alelFreq) <- unlist(lapply(auxped$markerdata,function(x){
                         res<-attr(x,"name");
                         return(res)
                         }))

 knownIds <- knownIds[!knownIds%in%QP]

 m <- c()
 for(i in seq_along(auxped$markerdat)){
  maux <-  apply(auxped$markerdata[[i]][knownIds,,drop=FALSE],1,function(x){
                    return(attr(auxped$markerdata[[i]],"alleles")[x])
                  })

  for(k in seq_along(maux)){
    if(length(maux[[k]])==0) maux[[k]]<-c(NA,NA)
  }
  tmaux <- (matrix(unlist(maux),ncol=2,byrow=TRUE))
  m <- cbind(m,tmaux)
 }
 colnames(m)<-paste(rep(systems,each=2),c("p","m"),sep="_")
 rownames(m)<-knownIds

 if(TRUE){#!exists("mmodel")){
  mmodel <- list()
  mmodel[["hombre"]]<-mmodel[["mujer"]]<-list(none=c())
 }else{
   warning("Check mut. model specification!\n")
   if(!is.list(mmodel)){
   if(!mmodel%in%MUT_MODEL_NAMES) warning(paste(mmodel,"not in", paste(MUT_MODEL_NAMES,collapse="|")))
   mmodel <- list(hombre=mmodel[1],mujer=mmodel[1])
  }
 }
 return(list(ped=auxped, markerEvidence=m,markerLinkage=markerLinkage,alelFreq=alelFreq,mmodel=mmodel))
}
