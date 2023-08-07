#' reverseSpit: a function for formatting.
#'
#' @param inList input for formatting.
#' @import paramlink
#' @import graphics
#' @export
#' @return A bayesian network.

reverseSplit<-function (inList) {
    if (length(inList) == 0) {
        return(inList)
    }
    lens = sapply(inList, length)
    nms = rep(names(inList), lens)
    vals = unlist(inList)
    split(nms, vals)
}
