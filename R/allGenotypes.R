#' Genotype combinations
#'
#' Auxiliary functions computing possible genotype combinations in a pedigree.
#' These are not normally intended for end users.
#'
#'
#' @param x a linkdat object.
#' @param partialmarker a marker object compatible with \code{x}.
#' @param ids a numeric with ID labels of one or more pedigree members.
#' @param chrom a character, either 'X' or 'AUTOSOMAL'. If missing, the 'chrom'
#'   attribute of \code{partialmarker} is used. If this is also missing, then
#'   'AUTOSOMAL' is taken as the default value.
#' @param make.grid a logical. If FALSE, a list is returned, otherwise
#'   \code{fast.grid} is applied to the list before returning it.
#' @param n a positive integer.
#' @param argslist a list of vectors.
#' @param as.list if TRUE, the output is a list, otherwise a matrix.
#' @return \code{allGenotypes} returns a matrix with 2 columns and \code{n +
#'   n*n(n-1)/2} rows containing all possible (unordered) genotypes at a
#'   biallelic locus with alleles \code{1,2,\dots{},n}. \code{fast.grid} is
#'   basically a stripped down version of expand.grid.
#'
#' @examples
#'
#' m = allGenotypes(2)
#' stopifnot(m == rbind(c(1,1), c(2,2), 1:2))
#'
#' @export
allGenotypes = function(n) rbind(cbind(seq_len(n), seq_len(n)), .comb2(n))

#' @rdname allGenotypes
#' @export
fast.grid = function(argslist, as.list = FALSE) {
    nargs = length(argslist)
    orep = nr = prod(lengths(argslist))
    if (nargs == 0L || nr == 0L)
        return(matrix(ncol = 0, nrow = 0))

    rep.fac = 1L
    res = NULL
    for (x in argslist) {
        nx = length(x)
        orep = orep/nx
        res = c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)])  #this is res[, i]
        rep.fac = rep.fac * nx
    }
    dim(res) = c(nr, nargs)
    if (as.list)
        res = lapply(seq_len(nr), function(r) res[r, ])
    res
}


#' @rdname allGenotypes
#' @export
geno.grid.subset = function(x, partialmarker, ids, chrom, make.grid = T) {
    int.ids = .internalID(x, ids)
    nall = attr(partialmarker, "nalleles")
    mutations = !is.null(attr(partialmarker, "mutmat"))
    if (missing(chrom))
        chrom = if (identical(attr(partialmarker, "chrom"), 23))
            "X" else "AUTOSOMAL"

    allg = allGenotypes(nall)
    allg_ref = 1000 * (allg[, 1] + allg[, 2]) + abs(allg[, 1] - allg[, 2])

    match_ref_rows = function(genomatr) {
        # In: matrix with 2 rows (each column a genotype). Out: vector of 'allg' row numbers
        sort.int(unique.default(match(1000 * (genomatr[1, ] + genomatr[2, ]) + abs(genomatr[1,
            ] - genomatr[2, ]), allg_ref)))
    }
    switch(chrom, AUTOSOMAL = {
        glist = .build_genolist(x, partialmarker, eliminate = ifelse(mutations, 0, 100))
        if (attr(glist, "impossible")) stop("Impossible partial marker")
        rows = lapply(glist[int.ids], match_ref_rows)
    }, X = {
        SEX = x$pedigree[, "SEX"]
        glist = .build_genolist_X(x, partialmarker, eliminate = ifelse(mutations, 0, 100))
        if (attr(glist, "impossible")) stop("Impossible partial marker")
        rows = lapply(int.ids, function(i) switch(SEX[i], glist[[i]], match_ref_rows(glist[[i]])))
    })
    if (make.grid)
        fast.grid(rows) else rows
}
.make.grid.subset = geno.grid.subset
