#' Marker functions
#'
#' Functions for setting and manipulating marker genotypes for 'linkdat'
#' objects compatible with fbnet. It was extracted from paramlink package, no longer maintained.
#'
#' @param x a \code{\link{linkdat}} object
#' @param ...  an even number of vectors, indicating individuals and their
#'   genotypes. See examples.
#' @param allelematrix a matrix with one row per pedigree member and two columns
#'   per marker, containing the alleles for a single marker.
#' @param m a \code{marker} object or a matrix with alleles. (In
#'   \code{SetMarkersfb} this matrix can contain data of several markers.)
#' @param missing a numeric - or character - of length 1, indicating the code
#'   for missing alleles.
#' @param chrom NA or an integer (the chromosome number of the marker).
#' @param pos NA or a non-negative real number indicating the genetic position
#'   (in cM) of the marker.
#' @param name NA or a character (the name of the marker).
#' @param mutmat a mutation matrix, or a list of two such matrices named
#'   'female' and 'male'. The matrix/matrices must be square, with the allele
#'   labels as dimnames, and each row must sum to 1 (after rounding to 3
#'   decimals).
#' @param annotations a list of marker annotations.
#' @param markernames NULL or a character vector.
#' @param chroms NULL or a numeric vector of chromosome numbers.
#' @param fromPos,toPos NULL or a single numeric.
#' @param marker,markers a numeric indicating which marker(s) to use/modify.
#' @param ids a numeric indicating individual(s) to be modified. In
#'   \code{swapGenotypesfb} this must have length 2.
#' @param genotype a vector of length 1 or 2, containing the genotype to be
#'   given the \code{ids} individuals. See examples.
#' @param alleles a numeric or character vector containing allele names.
#' @param afreq a numerical vector with allele frequencies. An error is given if
#'   they don't sum to 1 (rounded to 3 decimals).
#' @param new.alleles a numerical matrix of dimensions \code{length(ids),
#'   2*x$nMark}. Entries refer to the internal allele numbering.
#'
#' @return The \code{marker} function returns an object of class \code{marker}:
#'   This is a numerical 2-column matrix with one row per individual, and
#'   attributes 'alleles' (a character vector with allele names), 'nalleles'
#'   (the number of alleles) and 'missing' (the input symbol for missing marker
#'   alleles), 'chrom' (chromosome number), 'name' (marker identifier), 'pos'
#'   (chromosome position in cM).
#'
#'   For \code{addMarker}, \code{SetMarkersfb}, \code{removeMarkersfb},
#'   \code{modifyMarkerfb}, \code{modifyMarkerfbMatrix} and \code{swapGenotypesfb}, a
#'   \code{linkdat} object is returned, whose \code{markerdata} element has been
#'   set/modified.
#'
#'   For \code{getMarkersfb} a numeric vector containing marker numbers (i.e.
#'   their indices in \code{x$markerdata}) for the markers whose 'name'
#'   attribute is contained in \code{markernames}, 'chrom' attribute is
#'   contained in \code{chroms}, and 'pos' attribute is between \code{from} and
#'   \code{to}. NULL arguments are skipped, so \code{getMarkersfb(x)} will return
#'   \code{seq_len(x$nMark)} (i.e. all markers).
#'
#' @seealso \code{\link{linkdat}}
#'
#'
#'
#' @name markers
NULL

#' @rdname markers
#' @export
markerfb = function(x, ..., allelematrix, alleles = NULL, afreq = NULL, missing = 0, chrom = NA,
    pos = NA, name = NA, mutmat = NULL) {
    arglist = list(...)
    n = length(arglist)
    if (n == 0) {
        if (missing(allelematrix))
            m = matrix(missing, ncol = 2, nrow = x$nInd) else {
            m = as.matrix(allelematrix)
            stopifnot(nrow(m) == x$nInd, ncol(m) == 2)
        }
    }
    if (n > 0) {
        if (!missing(allelematrix))
            stop("Individual genotypes ar not allowed when 'allelematrix' is given.")
        if (n == 1 && length(arglist[[1]]) > 2)
            stop("Syntax error. See ?marker.")
        if (n > 1 && n%%2 != 0)
            stop("Wrong number of arguments.")

        fill = {
            if (n == 1)
                arglist[[1]] else missing
        }
        m = matrix(fill, ncol = 2, nrow = x$nInd, byrow = TRUE)  #create marker matrix with all individuals equal.

        for (i in (seq_len(n/2) * 2 - 1)) {
            # odd numbers up to n-1.
            ids = arglist[[i]]
            geno = arglist[[i + 1]]
            for (j in match(ids, x$orig.ids)) m[j, ] = geno
        }
    }

    .createMarkerObject(m, missing = missing, alleles = alleles, afreq = afreq, chrom = chrom,
        pos = pos, name = name, mutmat = mutmat)
}

#' @rdname markers
#' @export
addMarkerfb = function(x, m, ...) {
    if (is.matrix(m) || is.data.frame(m))
        stopifnot(nrow(m) == x$nInd, ncol(m) == 2)
    if (inherits(m, "marker"))
        m = list(m)
    if (is.list(m) && all(sapply(m, inherits, what = "marker")))
        return(SetMarkersfb(x, structure(c(x$markerdata, m), class = "markerdata")))
    if (!is.list(m) && length(m) == 1)
        m = matrix(m, ncol = 2, nrow = x$nInd)  #gives a nice way of setting an empty or everywhere-homozygous marker, e.g.: x=addMarker(x,0)
    mm = .createMarkerObject(m, ...)
    SetMarkersfb(x, structure(c(x$markerdata, list(mm)), class = "markerdata"))
}

#' @rdname markers
#' @export
SetMarkersfb = function(x, m, annotations = NULL, missing = 0) {
    if (is.null(m))
        markerdata_list = NULL
    else if (inherits(m, "marker"))
        markerdata_list = structure(list(m), class = "markerdata")
    else if (is.list(m) && all(sapply(m, inherits, what = "marker")))
        markerdata_list = structure(m, class = "markerdata")
    else if (inherits(m, "markerdata"))
        markerdata_list = m
    else if ((n <- ncol(m <- as.matrix(m))) == 0)
        markerdata_list = NULL
    else {
        if (is.character(m[1, 1]) && ("/" %in% strsplit(m[1, 1], "")[[1]])) {
            # if alleles are merged to 1 genotype per column
            splitvec = unlist(strsplit(m, "/", fixed = T))
            nrows = nrow(m)
            msplit = matrix(missing, ncol = 2 * n, nrow = nrows)
            msplit[, 2 * seq_len(n) - 1] = splitvec[2 * seq_len(n * nrows) - 1]
            msplit[, 2 * seq_len(n)] = splitvec[2 * seq_len(n * nrows)]
            m = msplit
        }
        if (ncol(m)%%2 != 0)
            stop("Uneven number of marker allele columns")
        nMark = ncol(m)/2
        if (!is.null(annotations)) {
            if (length(annotations) == 2 && !is.null(names(annotations)))
                annotations = rep(list(annotations), nMark)  # if given attrs for a single marker
            else if (length(annotations) != nMark)
                stop("Length of marker annotation list does not equal number of markers.")

            markerdata_list = lapply(1:nMark, function(i) {
                if (is.null(attribs <- annotations[[i]]))
                  return(NULL)
                mi = m[, c(2 * i - 1, 2 * i), drop = FALSE]
                do.call(.createMarkerObject, c(list(matr = mi, missing = missing), attribs))
            })
        } else {
            markerdata_list = lapply(1:nMark, function(i) {
                mi = m[, c(2 * i - 1, 2 * i), drop = FALSE]
                .createMarkerObject(mi, missing = missing)
            })
        }
        markerdata_list[sapply(markerdata_list, is.null)] = NULL
        class(markerdata_list) = "markerdata"
    }
    x$nMark = length(markerdata_list)
    x$markerdata = markerdata_list
    x$available = if (x$nMark > 0) x$orig.ids[rowSums(do.call(cbind, markerdata_list)) > 0] else numeric(0)
    x
}

.createMarkerObject = function(matr, name = NA, chrom = NA, pos = NA, alleles = NULL, afreq = NULL,
    mutmat = NULL, missing = 0) {
    if (is.null(alleles)) {
        vec = as.vector(matr)
        alleles = unique.default(vec[vec != missing])
        if (length(alleles) == 0)
            alleles = 1
    }
    if (!is.numeric(alleles) && !any(grepl("[^0-9\\.]", alleles)))
        alleles = as.numeric(alleles)
    all_ord = order(alleles)
    alleles = alleles[all_ord]
    nalleles = length(alleles)
    if (is.null(afreq))
        afreq = rep.int(1, nalleles)/nalleles else {
        if (length(afreq) != nalleles)
            stop("Number of alleles don't match length of frequency vector")
        if (round(sum(afreq), 2) != 1)
            warning(paste("Allele frequencies for marker", name, " do not sum to 1:", paste(afreq,
                collapse = ", ")))
        afreq = afreq[all_ord]
    }
    if (!is.null(mutmat)) {
        stopifnot(is.list(mutmat) || is.matrix(mutmat))
        # If single matrix given: make sex specific list
        if (is.matrix(mutmat)) {
            mutmat = .checkMutationMatrix(mutmat, alleles)
            mutmat = list(male = mutmat, female = mutmat)
        } else {
            stopifnot(length(mutmat) == 2, setequal(names(mutmat), c("female", "male")))
            mutmat$female = .checkMutationMatrix(mutmat$female, alleles, label = "female")
            mutmat$male = .checkMutationMatrix(mutmat$male, alleles, label = "male")
        }
    }
    m_obj = match(matr, alleles, nomatch = 0)
    attributes(m_obj) = list(dim = dim(matr), name = name, chrom = chrom, pos = pos, nalleles = nalleles,
        alleles = as.character(alleles), afreq = afreq, mutmat = mutmat, missing = missing,
        class = "marker")
    m_obj
}

.checkMutationMatrix = function(mutmat, alleles, label = "") {
    # Check that mutation matrix is compatible with allele number / names.  Sort matrix
    # according to the given allele order (this is important since dimnames are NOT used in
    # calculations).
    N = length(alleles)
    if (label != "")
        label = sprintf("%s ", label)
    if (any((dm <- dim(mutmat)) != N))
        stop(sprintf("Dimension of %smutation matrix (%d x %d) inconsistent with number of alleles (%d).",
            label, dm[1], dm[2], N))
    if (any(round(.rowSums(mutmat, N, N), 3) != 1))
        stop(sprintf("Row sums of %smutation matrix are not 1.", label))
    alleles.char = as.character(alleles)
    if (!setequal(rownames(mutmat), alleles) || !setequal(colnames(mutmat), alleles))
        stop(sprintf("Dimnames of %smutation do not match allele names.", label))
    m = mutmat[alleles.char, alleles.char]

    # lumbability: always lumpable (any alleles) if rows are identical (except diagonal)
    attr(m, "lumpability") = NA
    if (N > 2) {
        if (all(vapply(1:N, function(i) diff(range(m[-i, i])) == 0, logical(1))))
            attr(m, "lumpability") = "always"  # If each column has 0 range
    }
    m
}


.prettyMarkers = function(m, alleles = NULL, sep = "", missing = NULL, singleCol = FALSE, sex) {
    if (is.null(m))
        return(m)
    if (is.matrix(m))
        m = list(m)
    if ((n <- length(m)) == 0)
        return(m)
    if (is.null(alleles))
        alleles = lapply(m, attr, "alleles") else {
        if (!is.atomic(alleles))
            stop("The parameter 'alleles' must be NULL, or an atomic vector.")
        if (length(alleles) < max(unlist(lapply(m, attr, "nalleles"))))
            stop("The indicated 'alleles' vector has too few alleles for some markers.")
        alleles = rep(list(alleles), n)
    }
    if (is.null(missing))
        missing = unlist(lapply(m, attr, "missing")) else {
        if (!is.atomic(missing) || length(missing) != 1)
            stop("The parameter 'mising' must be NULL, or a numeric/character of length 1.")
        missing = rep(missing, n)
    }
    mNames = unlist(lapply(m, attr, "name"))
    mNames[is.na(mNames)] = ""
    pretty_m = do.call(c, lapply(seq_len(n), function(i) c(missing[i], alleles[[i]])[m[[i]] +
        1]))
    dim(pretty_m) = c(length(pretty_m)/(2 * n), 2 * n)
    if (singleCol) {
        al1 = pretty_m[, 2 * seq_len(n) - 1, drop = FALSE]
        al2 = pretty_m[, 2 * seq_len(n), drop = FALSE]
        m.matrix = matrix(paste(al1, al2, sep = sep), ncol = n)
        chrom_X = unlist(lapply(m, function(mm) identical(23L, as.integer(attr(mm, "chrom")))))
        if (any(chrom_X)) {
            males = (sex == 1)
            if (!all(hh <- al1[males, chrom_X] == al2[males, chrom_X]))
                warning("Male heterozygosity at X-linked marker detected.")
            m.matrix[males, chrom_X] = al1[males, chrom_X]
        }
        colnames(m.matrix) = mNames
        return(m.matrix)
    } else {
        nam = rep(mNames, each = 2)
        nam[nzchar(nam)] = paste(nam[nzchar(nam)], 1:2, sep = "_")
        colnames(pretty_m) = nam
        return(pretty_m)
    }
}

#' @rdname markers
#' @export
modifyMarkerfb = function(x, marker, ids, genotype, alleles, afreq, chrom, name, pos) {
    if (inherits(marker, "marker")) {
        if (nrow(marker) != x$nInd)
            stop("Wrong dimensions of marker matrix.")
        m = marker
    } else {
        if (!is.numeric(marker) || length(marker) != 1)
            stop("Argument 'marker' must be a single integer or an object of class 'marker'.")
        if (marker > x$nMark)
            stop("Indicated marker does not exist.")
        m = x$markerdata[[marker]]
    }
    mis = attr(m, "missing")

    if (!missing(alleles)) {
        stopifnot(is.atomic(alleles), is.numeric(alleles) || is.character(alleles))
        # if(is.numeric(alleles)) alleles = as.character(alleles)
        if (attr(m, "missing") %in% alleles)
            stop("The 'missing allele' character cannot be one of the alleles.")
        lena = length(alleles)
        if (lena == attr(m, "nalleles"))
            attr(m, "alleles") = as.character(alleles) else {
            num_als = unique.default(as.vector(m[m != 0]))
            if (lena < length(num_als))
                stop("Too few alleles.")
            pm = matrix(c(0, attr(m, "alleles"))[m + 1], ncol = 2)
            m = .createMarkerObject(pm, missing = 0, alleles = alleles, afreq = rep(1, lena)/lena,
                chrom = attr(m, "chrom"), pos = attr(m, "pos"), name = attr(m, "name"), mutmat = attr(m,
                  "mutmat"))
        }
    }
    if (!missing(afreq)) {
        if (round(sum(afreq), 2) != 1)
            stop("The allele frequencies don't sum to 1.")
        if (length(afreq) != attr(m, "nalleles"))
            stop("The length of allele frequency vector doesn't equal the number of alleles.")
        if (is.null(names(afreq)))
            attr(m, "afreq") = afreq else if (all(names(afreq) %in% attr(m, "alleles")))
            attr(m, "afreq") = afreq[attr(m, "alleles")] else stop("The names of the frequency vector don't match the allele names.")
    }

    changegeno = sum(!missing(ids), !missing(genotype))
    if (changegeno == 1)
        stop("The parameters 'ids' and 'genotype' must either both be NULL or both non-NULL.")
    if (changegeno == 2) {
        if (!is.atomic(genotype) || length(genotype) > 2)
            stop("The 'genotype' parameter must be a numeric or character vector of length 1 or 2.")
        pm = .prettyMarkers(list(m))
        for (i in .internalID(x, ids)) pm[i, ] = genotype
        if (!all(as.vector(pm) %in% c(attr(m, "alleles"), mis)))
            stop("Unknown allele(s). Please specify allele names using the 'alleles' argument.")
        attribs = attributes(m)
        m = match(pm, attr(m, "alleles"), nomatch = 0)
        attributes(m) = attribs
    }
    if (!missing(chrom))
        attr(m, "chrom") = chrom
    if (!missing(name))
        attr(m, "name") = name
    if (!missing(pos))
        attr(m, "pos") = pos
    if (inherits(marker, "marker"))
        return(m) else {
        x$markerdata[[marker]] = m
        x
    }
}


#' @rdname markers
#' @export
getMarkersfb = function(x, markernames = NULL, chroms = NULL, fromPos = NULL, toPos = NULL) {
    mnos = seq_len(x$nMark)
    if (!is.null(markernames)) {
        if (length(markernames) == 0)
            return(numeric(0))
        mnos = mnos[match(markernames, unlist(lapply(x$markerdata, function(m) attr(m, "name"))),
            nomatch = 0)]
    }
    if (!is.null(chroms)) {
        if (length(chroms) == 0)
            return(numeric(0))
        mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, "chrom"))) %in% chroms]
    }
    if (!is.null(fromPos))
        mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, "pos"))) >= fromPos]
    if (!is.null(toPos))
        mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, "pos"))) <= toPos]
    mnos
}

#' @rdname markers
#' @export

removeMarkersfb = function(x, markers = NULL, markernames = NULL, chroms = NULL, fromPos = NULL,
    toPos = NULL) {
    if (is.null(markers))
        markers = getMarkersfb(x, markernames, chroms, fromPos, toPos)
    if (is.null(markers) || length(markers) == 0)
        return(x)
    m = x$markerdata
    m[markers] = NULL
    SetMarkersfb(x, m)
}

#' @rdname markers
#' @export
swapGenotypesfb = function(x, ids) {
    assert_that(length(ids) == 2)
    ids = .internalID(x, ids)
    y = as.matrix(x)
    y[ids, -(1:6)] = y[ids[2:1], -(1:6)]
    restore_linkdat(y)
}

#TODO!
#print.marker = function(m) {
#    cat(sprintf("Marker name: %s\n", attr(m, 'name')))
#    chrom = attr(m, 'chrom')
#    pos = attr(m, 'pos')
#    if(!is.na(chrom)) pos = paste(chrom, pos, sep=" - ")
#    cat(sprintf("Position: %s\n", pos))
#    cat("Alleles and frequencies:\n")
#    afreq = attr(m, 'afreq')
#    names(afreq) = attr(m, 'allele')
#    print(afreq)
#}

#' @rdname markers
#' @export
modifyMarkerfbMatrix = function(x, ids, new.alleles) {
    ids = .internalID(x, ids)
    y = as.matrix(x)
    y[ids, -(1:6)] = new.alleles
    restore_linkdat(y)
}

.setSNPfreqs = function(x, newfreqs) {
    stopifnot(all(vapply(x$markerdata, function(m) attr(m, "nalleles"), numeric(1)) == 2))
    newfreqs = rep(newfreqs, length = x$nMark)
    for (i in seq_len(x$nMark)) attr(x$markerdata[[i]], "afreq") = c(newfreqs[i], 1 - newfreqs[i])
    x
}
