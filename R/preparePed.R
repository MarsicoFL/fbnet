#' preparePed: a function for simulating genetic data from untyped individuals conditioned on known genotypes.
#'
#' @param ped A ped object with information of the genotyped members. The ped object must be in Familias format.
#' @param available Genotyped individuals IDs. 
#' @param lLociFreq Allele frequencies.
#' @param rseed Seed used for simulations.
#' @import paramlink
#' @import igraph
#' @export
#' @return A ped object.

preparePed<-function(ped,available,lLociFreq,rseed=NULL){

 if(!is.null(rseed)) set.seed(rseed)
 for(i in seq_along(lLociFreq)){
    alleles <- as.numeric(lLociFreq[[i]][,"Alelo"])
    afreq<-lLociFreq[[i]][,"freq"]/sum(lLociFreq[[i]][,"freq"])
    iout <- which(alleles==0)
    if(length(iout)>0){
     alleles<-alleles[-iout]
     afreq<-afreq[-iout]
    }

    paux   <- markerSim(ped, N=1, alleles=alleles,afreq=afreq,verbose=FALSE)
    mallele<- matrix(attr(paux$markerdata[[1]],"alleles")[paux$markerdata[[1]]],ncol=2)

    maux <- marker(ped,allelematrix=mallele,alleles=alleles,afreq=afreq,name=names(lLociFreq)[i])
    ped  <- addMarker(ped,maux)
 }

 inotEv <- ped$orig.ids[!ped$orig.ids%in%available]
 ped    <-removeEvidenceFromPed(ped,inotEv)
 return(ped)
}

