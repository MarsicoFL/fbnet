#' Initialize a Bayesian network from a pedtools pedigree
#'
#' Builds the list structure returned by \code{\link{initBN}} directly
#' from a \code{pedtools::ped} object, without routing through the
#' archived \code{paramlink} package. This is the recommended entry
#' point for modern workflows that construct pedigrees via
#' \code{pedtools}, \code{forrel}, or \code{pedsuite} helpers such as
#' \code{linearPed}, \code{cousinPed}, \code{avuncularPed}, or
#' \code{nuclearPed}.
#'
#' @param x A \code{pedtools::ped} object. Must have markers attached
#'   (via \code{pedtools::setMarkers}) with allele frequencies defined.
#' @param available Character or integer vector of IDs that are
#'   genotyped (have non-missing evidence). If \code{NULL} (the
#'   default), the set of typed individuals is inferred automatically
#'   from the pedigree: an individual is considered available if they
#'   have at least one non-missing genotype across the attached markers.
#' @param QP Character or integer vector of Query Person IDs (the
#'   individuals whose genotype posterior we want to query). Used by
#'   \code{\link{buildBN}} downstream. Can be \code{NULL} and supplied
#'   later.
#'
#' @return A list with the same shape as the return value of
#'   \code{\link{initBN}}, namely
#'   \describe{
#'     \item{\code{ped}}{An auxiliary pedigree object with fields
#'       \code{pedigree}, \code{ped}, \code{orig.ids}, \code{markerdata},
#'       \code{available}, and \code{nMark}, mimicking the subset of
#'       the \code{paramlink} linkdat interface that downstream fbnet
#'       functions consume.}
#'     \item{\code{markerEvidence}}{A character matrix with rows
#'       indexed by typed individuals and columns \code{<marker>_p} /
#'       \code{<marker>_m} giving the observed allele labels.}
#'     \item{\code{markerLinkage}}{A data frame with columns
#'       \code{linkedTo} and \code{recomb}; defaults to unlinked.}
#'     \item{\code{alelFreq}}{A named list of allele frequency vectors,
#'       one per marker.}
#'     \item{\code{mmodel}}{A nested list with the default no-mutation
#'       model for each sex.}
#'   }
#'
#' @details
#' The \code{paramlink} package was archived from CRAN in 2022, and its
#' replacement \code{pedtools} uses a different in-memory representation.
#' The original \code{\link{initBN.fromPed}} function requires a
#' \code{paramlink} linkdat object, which makes it difficult to use
#' fbnet with modern pedigree workflows. \code{from_pedtools} constructs
#' the equivalent list structure by hand from a \code{pedtools::ped}
#' object, so that the rest of the fbnet pipeline (\code{\link{buildBN}},
#' \code{\link{buildCPTs}}, \code{\link{imposeEvidence}}, etc.) can be
#' run on pedigrees built with modern tools such as \code{forrel} and
#' \code{pedsuite}.
#'
#' @note
#' This function handles the \emph{entry point}. Some downstream fbnet
#' functions still use \code{@import paramlink} in their roxygen
#' declarations, so a fully \code{paramlink}-free fbnet pipeline may
#' still require additional refactoring. For basic use cases (single
#' marker, small pedigree), the return value of \code{from_pedtools}
#' can be passed directly to \code{\link{buildBN}}.
#'
#' @examples
#' \dontrun{
#'   library(pedtools)
#'   library(forrel)
#'
#'   # Build a 2nd-degree grandparent-grandchild pedigree with one marker
#'   ped <- linearPed(2)
#'   m <- marker(ped, afreq = c(A = 0.3, B = 0.4, C = 0.3), name = "M1")
#'   ped <- setMarkers(ped, m)
#'   ped <- setGenotype(ped, id = "1", marker = 1, geno = "A/B")
#'
#'   pbn <- from_pedtools(ped)
#'   str(pbn)
#' }
#'
#' @seealso \code{\link{initBN}}, \code{\link{initBN.fromPed}},
#'   \code{\link{buildBN}}.
#'
#' @export
from_pedtools <- function(x, available = NULL, QP = NULL) {
  if (!requireNamespace("pedtools", quietly = TRUE))
    stop("The 'pedtools' package is required for from_pedtools(). ",
         "Install it with install.packages('pedtools').")

  if (!inherits(x, "ped"))
    stop("`x` must be a pedtools::ped object (use pedtools::linearPed, ",
         "cousinPed, avuncularPed, or similar).")

  ids <- as.character(labels(x))
  n   <- length(ids)

  # Resolve parent IDs from FIDX/MIDX (index into the ID vector).
  fidx <- as.integer(x$FIDX)
  midx <- as.integer(x$MIDX)
  fid  <- ifelse(fidx == 0, "0", ids[pmax(fidx, 1)])
  mid  <- ifelse(midx == 0, "0", ids[pmax(midx, 1)])
  fid[fidx == 0] <- "0"
  mid[midx == 0] <- "0"

  sex_num <- as.integer(x$SEX)   # 1 = male, 2 = female in pedtools
  aff     <- rep("1", n)

  # Infer typed individuals if `available` is not provided
  if (is.null(available)) {
    available <- character(0)
    if (pedtools::nMarkers(x) > 0) {
      typed <- vapply(ids, function(id) {
        for (k in seq_len(pedtools::nMarkers(x))) {
          g <- pedtools::genotype(x, id = id, marker = k)
          if (!is.null(g) && length(g) > 0 &&
              !all(is.na(g)) &&
              !all(g %in% c("-", "0", "", NA_character_))) return(TRUE)
        }
        FALSE
      }, logical(1))
      available <- ids[typed]
    }
  } else {
    available <- as.character(available)
  }

  # Pedigree matrix in paramlink linkdat shape
  pped_mat <- cbind(
    ID  = ids,
    FID = fid,
    MID = mid,
    SEX = as.character(sex_num),
    AFF = aff
  )

  # Build markerdata in a paramlink-compatible shape: each marker is an
  # integer matrix (n x 2) with attributes carrying alleles / freqs / name.
  markerdata <- list()
  alelFreq   <- list()
  for (k in seq_len(pedtools::nMarkers(x))) {
    m_k       <- x$MARKERS[[k]]
    mname     <- attr(m_k, "name") %||% paste0("M", k)
    mk_alleles <- attr(m_k, "alleles")
    mk_afreq   <- as.numeric(attr(m_k, "afreq"))
    names(mk_afreq) <- mk_alleles

    # pedtools markers are integer matrices n x 2 with indices into alleles.
    mat <- matrix(as.integer(m_k), nrow = n, ncol = 2)

    attr(mat, "name")     <- mname
    attr(mat, "chrom")    <- attr(m_k, "chrom") %||% NA
    attr(mat, "pos")      <- attr(m_k, "posMb") %||% NA
    attr(mat, "nalleles") <- length(mk_alleles)
    attr(mat, "alleles")  <- mk_alleles
    attr(mat, "afreq")    <- mk_afreq
    attr(mat, "missing")  <- 0L
    attr(mat, "mutmat")   <- NULL
    class(mat)            <- "marker"

    markerdata[[k]] <- mat
    alelFreq[[mname]] <- mk_afreq
  }

  systems <- names(alelFreq)
  nsys    <- length(systems)

  # Build the markerEvidence matrix: rows = typed individuals,
  # cols = <marker>_p, <marker>_m with observed allele labels.
  if (length(available) > 0 && nsys > 0) {
    m <- matrix(NA_character_, nrow = length(available), ncol = 2 * nsys)
    rownames(m) <- available
    colnames(m) <- as.vector(t(outer(systems, c("p", "m"), paste, sep = "_")))

    for (k in seq_len(nsys)) {
      marker_k <- markerdata[[k]]
      alleles_k <- attr(marker_k, "alleles")
      for (ai in seq_along(available)) {
        row_idx <- match(available[ai], ids)
        idx_p <- marker_k[row_idx, 1]
        idx_m <- marker_k[row_idx, 2]
        if (!is.na(idx_p) && idx_p > 0)
          m[ai, 2 * k - 1] <- alleles_k[idx_p]
        if (!is.na(idx_m) && idx_m > 0)
          m[ai, 2 * k] <- alleles_k[idx_m]
      }
    }
  } else {
    m <- matrix(character(0), nrow = 0, ncol = 0)
  }

  markerLinkage <- data.frame(
    linkedTo = rep("", nsys),
    recomb   = rep(0.5, nsys),
    stringsAsFactors = FALSE
  )
  rownames(markerLinkage) <- systems

  mmodel <- list()
  mmodel[["hombre"]] <- list(none = c())
  mmodel[["mujer"]]  <- list(none = c())

  # Build the auxped object with the minimal shape downstream fbnet
  # functions consume.
  auxped <- list(
    pedigree = pped_mat,
    ped      = pped_mat,
    orig.ids = ids,
    markerdata = markerdata,
    available  = available,
    nMark      = nsys
  )

  if (!is.null(QP)) {
    QP <- as.character(QP)
    qp_rows <- auxped$ped[, "ID"] %in% QP
    auxped$ped[qp_rows, "AFF"] <- "2"
    auxped$pedigree[qp_rows, "AFF"] <- "2"
  }

  list(
    ped            = auxped,
    markerEvidence = m,
    markerLinkage  = markerLinkage,
    alelFreq       = alelFreq,
    mmodel         = mmodel
  )
}

# Helper operator for `%||%` (fallback for missing attrs)
`%||%` <- function(a, b) if (is.null(a)) b else a
