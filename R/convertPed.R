#' Convert pedigree format from pedtools to paramlink
#'
#' Converts the `ped` object to a [paramlink::linkdat] object.
#'
#' @param x A `ped` object
#' @param verbose A logical
#'
#' @return A linkdat object:
#' @import pedtools
#'
#' @export
convertPed = function(x, verbose=F) {
  if (!requireNamespace("paramlink", quietly = TRUE))
    stop2("Package 'paramlink' is not installed")

  # famid = x$FAMID
  # if(famid == "") famid = 1
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
