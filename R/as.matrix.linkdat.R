#' linkdat to matrix conversion
#'
#' Converts a \code{linkdat} object to a matrix (basically following a
#' pre-makeped LINKAGE format), with marker annotations and other info attached
#' as attributes.
#'
#' \code{restore_linkdat} is the reverse of \code{as.matrix}.
#'
#' The way \code{linkdat} objects are created in paramlink, marker data
#' are stored as a list of \code{marker} objects. Each of these is essentially a
#' matrix with various attributes like allele frequencies, map info a.s.o.. This
#' format works well for marker-by-marker operations (e.g. likelihoods and LOD
#' scores), but makes it somewhat awkward to operate 'horizontally', i.e.
#' individual-by-individual, for instance if one wants to delete all genotypes
#' of a certain individual, or rearrange the pedigree in some way.
#'
#' It is therefore recommended to convert the \code{linkdat} object to a matrix
#' first, do the necessary manipulations on the matrix, and finally use
#' \code{restore_linkdat}. Attributes are often deleted during matrix
#' manipulation, so it may be necessary to store them in a variable and feed
#' them manually to \code{restore_linkdat} using the \code{attrs} argument.
#'
#' With default parameters, \code{restore_linkdat(as.matrix(x))} should
#' reproduce \code{x} exactly.
#'
#' @param x a \code{linkdat} object. In \code{restore_linkdat}: A
#'   numerical matrix in LINKAGE format.
#' @param include.attrs a logical indicating if marker annotations and other
#'   info should be attached as attributes. See value.
#' @param attrs a list containing marker annotations and other \code{linkdat}
#'   info compatible with \code{x}, in the format produced by \code{as.matrix}.
#'   If NULL, the attributes of \code{x} itself are used.
#' @param checkped a logical, forwarded to \code{linkdat}. If FALSE, no
#'   checks for pedigree errors are performed.
#' @param \dots not used.
#'
#' @return For \code{as.matrix}: A matrix with \code{x$nInd} rows and \code{6 +
#'   2*x$nMark} columns.  The 6 first columns describe the pedigree in LINKAGE
#'   format, and the remaining columns contain marker alleles, using the
#'   internal (numerical) allele coding and 0 for missing alleles. If
#'   \code{include.attrs = TRUE} the matrix has the following attributes:
#'   \itemize{ \item \code{markerattr} (a list of marker annotations) \item
#'   \code{available} (the availability vector) \item \code{model} (the disease
#'   model, if present) \item \code{plot.labels} (plot labels, if present) \item
#'   \code{orig.ids} (original individual IDs) }
#'
#'   For \code{restore_linkdat}: A \code{linkdat} object.
#'
#' @seealso \code{linkdat}, \code{as.data.frame.linkdat}
#'
#'
#' @export
as.matrix.linkdat = function(x, include.attrs = TRUE, ...) {
    p = do.call(cbind, c(list(FAMID = x$famid, relabelfb(x$pedigree, x$orig.ids)), x$markerdata))
    if (include.attrs) {
        attr(p, "markerattr") = lapply(x$markerdata, attributes)
        attr(p, "available") = x$available
        attr(p, "model") = x$model
        attr(p, "plot.labels") = x$plot.labels
        attr(p, "orig.ids") = x$orig.ids  # needed to restore plot.labels
    }
    p
}

#' @rdname as.matrix.linkdat
#' @export
restore_linkdat = function(x, attrs = NULL, checkped = TRUE) {
    if (is.null(attrs))
        attrs = attributes(x)
    y = linkdat(x[, 1:6, drop = F], model = attrs$model, checkped = checkped, verbose = FALSE)

    # Create marker objects
    markers = x[, -(1:6), drop = F]
    nMark = ncol(markers)/2
    if (nMark == 0)
        markerdata_list = NULL
    else {
        markerattr = attrs$markerattr
        markerdata_list = lapply(seq_len(nMark), function(k) {
            m = markers[, c(2 * k - 1, 2 * k), drop = F]
            attributes(m) = c(markerattr[[k]][-1], list(dim = dim(m)))
            m
        })
        class(markerdata_list) = "markerdata"
    }

    y = SetMarkersfb(y, markerdata_list)
    y = setAvailable(y, intersect(attrs$available, y$orig.ids))
    if (!is.null(pl <- attrs$plot.labels)) {
        y$plot.labels = pl[match(y$orig.ids, attrs$orig.ids)]
        y$plot.labels[is.na(y$plot.labels)] = ""
    }
    y
}

