#' runFamilias: a function for running Familias likelihood ratio computations.
#'
#' @param ped A ped object with information of the genotyped members. The ped object must be in Familias format.
#' @import Familias
#' @import paramlink
#' @import graphics
#' @export
#' @return Computation results.

runFamilias<-function(ped=NULL){
geno <- QP <- lLociFreq <- NULL
  if(!is.null(ped)){
    persons <- ped$orig.ids
    pid     <- ped$pedigree[,"FID"];pid[pid==0]<-NA
    mid     <- ped$pedigree[,"MID"];mid[mid==0]<-NA
    sex     <- ped$pedigree[,"SEX"];sex <- ifelse(sex==1,"male","female")
    systems <- unlist(lapply(ped$markerdata,function(x){attr(x,"name")}))
    bSimuData <- FALSE
    bplotped  <- FALSE
    datamatrix<-c()
    for(i in seq_along(systems)){
      m <- c(matrix(ped$markerdata[[i]],ncol=2))
      m[m!=0]<-attr(ped$markerdata[[i]],"alleles")[m[m!=0]]
      mallele<- matrix(m,ncol=2)
      colnames(mallele)<-paste(paste0("locus",i),c("p","m"),sep=".")
      rownames(mallele)<-ped$orig.ids
      datamatrix<-cbind(datamatrix,mallele)
      datamatrix[datamatrix=="0"]<-NA
    }
    datamatrix[QP,]<-t(apply(geno[QP,,drop=FALSE],1,function(x){unlist(strsplit(x,"/"))}))
  }

  ped1 <- FamiliasPedigree(id=persons,dadid=pid,momid=mid,sex=sex)
  if(!is.null(QP)){
    ppid <- pid
    mmid <- mid
    ppid[persons%in%QP] <- mmid[persons%in%QP] <- NA
    ped0 <- FamiliasPedigree(id=persons,dadid=ppid,momid=mmid,sex=sex)
  }

  myloci <- list()
  for(i in seq_along(systems)){
    if(!systems[i]%in%names(lLociFreq)) next
    freqs <- lLociFreq[[systems[i]]][,"freq"]
    anames<- lLociFreq[[systems[i]]][,"Alelo"]
    anames[anames=="0"]<-"100"
    if(sum(freqs)>1) freqs<-freqs/sum(freqs)
    if((1-sum(freqs))> 1e-10){
      freqs<-c(freqs,1-freqs)
      anames<-c(anames,"ExtraAlelle")
    }
    locus<-FamiliasLocus(frequencies=freqs,
                         allelenames=anames,
                         name=systems[i])
    myloci[[systems[i]]]<-locus
  }


  auxped<-Familias2linkdat(ped1,NULL,myloci)
  if(bSimuData){
    knownIds <- c("1","3","4","6")
    datamatrix<-c()
    for(i in seq_along(systems)){
      afreq <- myloci[[i]]$alleles

      paux<-markerSim(auxped,N=1,available=knownIds,alleles=names(afreq),afreq=afreq,seed=123457,verbose=FALSE)
      m <- c(matrix(paux$markerdata[[1]],ncol=2))
      m[m!=0]<-attr(paux$markerdata[[1]],"alleles")[m[m!=0]]
      mallele<- matrix(m,ncol=2)
      colnames(mallele)<-paste(paste0("locus",i),c("p","m"),sep=".")
      rownames(mallele)<-paux$orig.ids
      maux   <- marker(auxped,allelematrix=mallele,alleles=names(afreq),afreq=afreq,name=systems[i])
      auxped <- addMarker(auxped,maux)

      datamatrix<-cbind(datamatrix,mallele)
      datamatrix[datamatrix=="0"]<-NA

    }
  }else{

    rownames(datamatrix)<-persons
    knownIds <- rownames(datamatrix)[!is.na(datamatrix[,1])]
    for(i in seq_along(systems)){
      if(!systems[i]%in%names(lLociFreq)) next
      freqs <- lLociFreq[[systems[i]]][,"freq"]
      anames<- lLociFreq[[systems[i]]][,"Alelo"]
      anames[anames=="0"]<-"100"
      if(sum(freqs)>1) freqs<-freqs/sum(freqs)
      if((1-sum(freqs))>1e-4){
        freqs<-c(freqs,1-freqs)
        anames<-c(anames,"ExtraAlelle")
      }

      maux   <- marker(auxped,allelematrix=datamatrix[,(1:2)+2*(i-1)],alleles=anames,afreq=freqs,name=systems[i])
      auxped <- addMarker(auxped,maux)
    }

  }
  if(!is.null(QP)){
    auxped$pedigree[auxped$orig.ids%in%QP,"AFF"]<-2
  }
  if(bplotped)  plot(auxped,mark=1:min(3,length(myloci)))

  if(FALSE){
    auxped <- Familias2linkdat(ped1,NULL,myloci)
    for(i in seq_along(systems)){
      freqs <- lLociFreq[[systems[i]]][,"freq"]
      anames<- lLociFreq[[systems[i]]][,"Alelo"]
      anames[anames=="0"]<-"100"
      if(sum(freqs)>1) freqs<-freqs/sum(freqs)
      if(sum(freqs)<1){
        freqs<-c(freqs,1-freqs)
        anames<-c(anames,"ExtraAlelle")
      }

      maux   <- marker(auxped,allelematrix=datamatrix[,(1:2)+2*(i-1)],alleles=anames,afreq=freqs,name=systems[i])
      auxped <- addMarker(auxped,maux)
    }


    ppid<-pid;ppid[is.na(ppid)]<-0
    mmid<-mid;mmid[is.na(mmid)]<-0
    ddatamatrix<-datamatrix
    ddatamatrix[is.na(ddatamatrix)]<-0
    df <- data.frame(ID=persons,FID=ppid,MID=mmid,SEX=ifelse(sex=="male",1,2),AFF=rep(1,length(mmid)),ddatamatrix)


    ld1<-linkdat(df,annotation=mmyloci)
  }

  if(!is.null(QP)){

    saux <- paste(rep(names(myloci),each=2),c("p","m"),sep=".")
    icol <- which(colnames(datamatrix)%in%saux)
    if(length(icol)==0){
      saux <- paste(rep(paste0("locus",seq_along(myloci)),each=2),c("p","m"),sep=".")
      icol <- which(colnames(datamatrix)%in%saux)
    }

    L1 <- FamiliasPosterior(ped1,myloci,datamatrix[,icol])$likelihoodsPerSystem

    datamatrix0 <- datamatrix[,icol]
    datamatrix0[QP,]<-c(NA,NA)

    lrmp<-list()
    for(i in seq_along(myloci)){
      locus<-names(myloci)[i]
      pp <- 1
      for(j in seq_along(QP)){
        al <- as.character(unlist(datamatrix[QP[j],,drop=T]))
        al <- al[((2*(i-1))+1):(2*i)]
        pp <- pp * prod(myloci[[i]]$alleles[al])*2^(al[1]!=al[2])
      }
      lrmp[[myloci[[i]]$locusname]] <- pp
    }

    L0 <- FamiliasPosterior(ped1,myloci,datamatrix0)$likelihoodsPerSystem
    L0 <- L0 * unlist(lrmp[rownames(L0)])

    res <- cbind(L0=L0,L1=L1,LR=L1/L0)
    colnames(res)<-c("L0","L1","LR")
  }else{
    res <- FamiliasPosterior(list(ped1),myloci,datamatrix)
  }

  return(res)

}
