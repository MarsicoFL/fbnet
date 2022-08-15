#' initBN: a function to initialize the bayesian network.
#'
#' @param ped A ped object with information of the genotyped members. The ped object must be in Familias format.
#' @param bplotped An alternative ped object to be compared. 
#' @import paramlink
#' @import graphics
#' @examples 
#' pbn  <- initBN(toyped)
#' @export
#' @return A bayesian network.

initBN<-function(ped=NULL,bplotped=FALSE){

MUT_MODEL_NAMES <- c("none","equal","stepwise")
 
#if(exists("mmodel")) rm(mmodel)

 if(is.null(ped)){
   res <- initBN.fromVars(bplotped=bplotped)
 }else{
   res <- initBN.fromPed(ped,bplotped=bplotped)
 }
 return(res)
}

