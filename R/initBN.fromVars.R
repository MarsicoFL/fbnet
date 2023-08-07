#' initBN.fromVars: a function to initialize the bayesian network.
#'
#' @param bplotped An alternative ped object to be compared. 
#' @import paramlink
#' @import graphics
#' @export
#' @return A bayesian network.

initBN.fromVars<-function(bplotped){
ped1 <- FamiliasPedigree(id=persons,dadid=pid,momid=mid,sex=sex)

persons <- pid <- mid <- sex <- systems <- lLociFreq <- linkageR <- bSimuData <- simuKnownIds <- QP <- NULL
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

 
auxped<-paramlink::Familias2linkdat(ped1,NULL,myloci) 

if(bSimuData){ 
  knownIds <- simuKnownIds
  datamatrix<-c()
  for(i in seq_along(systems)){
   afreq <- myloci[[i]]$alleles
 
   paux<-paramlink::markerSim(auxped,N=1,available=knownIds,alleles=names(afreq),afreq=afreq,seed=123457,verbose=FALSE)
   m <- c(matrix(paux$markerdata[[1]],ncol=2))
   m[m!=0]<-attr(paux$markerdata[[1]],"alleles")[m[m!=0]]
   mallele<- matrix(m,ncol=2)
   colnames(mallele)<-paste(paste0("locus",i),c("p","m"),sep=".")
   rownames(mallele)<-paux$orig.ids
   maux   <- paramlink::marker(auxped,allelematrix=mallele,alleles=names(afreq),afreq=afreq,name=systems[i])
   auxped <- paramlink::addMarker(auxped,maux)
   
   datamatrix<-cbind(datamatrix,mallele)
   datamatrix[datamatrix=="0"]<-NA   
   
  }    
 }else{

  rownames(datamatrix)<-persons              
  knownIds <- rownames(datamatrix)[!is.na(datamatrix[,1])]
  for(i in seq_along(systems)){
   freqs <- lLociFreq[[systems[i]]][,"freq"]
   anames<- lLociFreq[[systems[i]]][,"Alelo"]
   anames[anames=="0"]<-"100"
   if(sum(freqs)>1) freqs<-freqs/sum(freqs)
   if(1-sum(freqs)>1e-4){ 
    freqs<-c(freqs,1-freqs)
    anames<-c(anames,"ExtraAlelle")
   }   
  
   maux   <- paramlink::marker(auxped,allelematrix=datamatrix[,(1:2)+2*(i-1)],alleles=anames,afreq=freqs,name=systems[i])
   auxped <- paramlink::addMarker(auxped,maux)
  } 
  
} 
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
                         
 m <- datamatrix[knownIds,]
 colnames(m)<-paste(rep(systems,each=2),c("p","m"),sep="_")
 markerEvidence<-t(apply(m,1,as.character))
 colnames(markerEvidence)<-colnames(m)
 rownames(markerEvidence)<-rownames(m)
 
 mmodel <- list()
 mmodel[["hombre"]]<-list(none=c())
 mmodel[["mujer"]]<-list(none=c())
 
  return(list(ped=auxped, markerEvidence=markerEvidence,markerLinkage=markerLinkage,alelFreq=alelFreq,
              mmodel=mmodel))
 
} 

