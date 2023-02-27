#' getConditional: a function for obtaining the coditional probability tables based on a given evidence.
#'
#' @param lf A list of joint probabilities.
#' @export
#' @return A list of conditioned probabilities.

getConditional<-function(lf){
     lf<-lapply(lf,function(x){
                    x$prob <- x$prob/sum(x$prob)
                    return(x)
                  })
      return(lf)
 }
