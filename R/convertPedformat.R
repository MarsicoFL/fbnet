#' convertPedformat: a function for converting a pedtools ped onject to a famlink ped object.
#'
#' @param x A pedtools ped object.
#' @param verbose  Function output.
#' @import paramlink
#' @export
#' @return A dataframe with LRs.

convertPedformat = function(x, verbose=FALSE) {
  afreq <- alleles <- chrom <- isXmarker <- mutmod <- nAlleles <- name <- posMb <- NULL
  famid = 1

  mlist = x$MARKERS

  x$MARKERS = NULL
  p = cbind(famid, as.matrix(x), 1)
  colnames(p) = c("FAMID", "ID", "FID", "MID", "SEX", "AFF")

  y = paramlink::linkdat(p, verbose=verbose)

  if(!is.null(mlist)) {
    mlist = lapply(mlist, function(m) {
      attributes(m) =
        list(dim = dim(m),
             name = name(m),
             chrom = if(isXmarker(m)) 23 else chrom(m),
             pos = posMb(m),
             nalleles = nAlleles(m),
             alleles = alleles(m),
             afreq = as.vector(afreq(m)),
             missing = 0,
             mutmat = mutmod(m),
             class = "marker")
      m
    })
    y = paramlink::setMarkers(y, mlist)
  }
  y
}
