#' Is an object a linkdat object?
#'
#' Functions for checking whether an object is a \code{\link{linkdat}} object, a
#' \code{\link{singletonfb}} or a list of such.
#'
#' Note that the \code{singletonfb} class inherits from \code{linkdat}, so if
#' \code{x} is a singletonfb, \code{is.linkdat(x)} returns TRUE.
#'
#' @param x Any \code{R} object.
#' @return For \code{is.linkdat}: TRUE if \code{x} is a linkdat (or singletonfb)
#'   object, and FALSE otherwise.\cr For \code{is.singletonfbfb}: TRUE if \code{x}
#'   is a singletonfb object, and FALSE otherwise.\cr For \code{is.linkdat.list}:
#'   TRUE if \code{x} is a list of linkdat/singletonfb objects.
#' @seealso \code{\link{linkdat}}
#'
#' @export
is.linkdat = function(x)
    inherits(x, "linkdat")

#' @rdname is.linkdat
#' @export
is.singletonfbfb = function(x)
    inherits(x, "singletonfb")

#' @rdname is.linkdat
#' @export
is.linkdat.list = function(x)
    isTRUE(is.list(x) && all(sapply(x, inherits, "linkdat")))

