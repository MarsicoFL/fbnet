#' Convert 'Familias' output to linkdat objects
#' 
#' This function is extracted from paramlik package, not currently maintained.
#' Familias is a widely used program for computations in forensic genetics. The
#' function documented here facilitates the use of \code{paramlink} for
#' specialized computations which are not implemented in Familias, e.g.
#' conditional simulations.
#'
#' The Familias program represents pedigrees and marker data in a way that
#' differs from \code{paramlink} in several ways, mostly because of
#' \code{paramlink}'s stricter definition of a 'pedigree'. In \code{paramlink},
#' a pedigree must be connected, have numerical IDs, and each member must have
#' either 0 or 2 parentsfb present in the pedigree. None of this is required by
#' \code{FamiliasPedigree} objects. The conversion function
#' \code{Familias2linkdat} takes care of all of these potential differences: It
#' converts each \code{FamiliasPedigree} into a list of connected linkdat objects,
#' additional parentsfb are added where needed, and non-numerical ID labels are
#' stored in the \code{plot.labels} slot of the linkdat object(s).
#'
#' @param familiasped A \code{FamiliasPedigree} object or a
#'   list of such.
#' @param datamatrix A data frame with two columns per marker (one for each
#'   allele) and one row per individual.
#' @param loci A \code{FamiliasLocus} object or a list of such.
#' @param ID An integer vector: Individual ID.
#' @param FID An integer vector: ID of father.
#' @param MID An integer vector: ID of mother.
#' @return A \code{\link{linkdat}} object, or a list of such.
#' @author Magnus Dehli Vigeland, Thore Egeland
#' @references Windows Familias is freely available from
#'   \url{https://familias.name}.
#' @export
Familias2linkdat = function(familiasped, datamatrix, loci) {

    ### If first argument is a list of FamiliasPedigrees, convert one at a time.
    if (is.list(familiasped) && inherits(familiasped[[1]], "FamiliasPedigree")) {
        res = lapply(familiasped, function(p) Familias2linkdat(p, datamatrix = datamatrix,
            loci = loci))
        return(res)
    }

    ### Part 1: pedigree

    K = length(familiasped$id)
    ped = data.frame(ID = 1:K, FID = familiasped$findex, MID = familiasped$mindex, SEX = ifelse(familiasped$sex ==
        "male", 1, 2), AFF = 1)

    fatherMissing = which(ped$FID == 0 & ped$MID != 0)
    motherMissing = which(ped$FID != 0 & ped$MID == 0)

    nFath = length(fatherMissing)
    nMoth = length(motherMissing)

    newFathers = seq(K + 1, length = nFath)
    newMothers = seq(K + 1 + nFath, length = nMoth)

    plotLabels = c(familiasped$id, paste("added", seq_len(nFath + nMoth)))

    # add new fathers
    if (nFath > 0) {
        ped = rbind(ped, data.frame(ID = newFathers, FID = 0, MID = 0, SEX = 1, AFF = 1))
        ped[fatherMissing, "FID"] = newFathers
    }

    # add new mothers
    if (nMoth > 0) {
        ped = rbind(ped, data.frame(ID = newMothers, FID = 0, MID = 0, SEX = 2, AFF = 1))
        ped[motherMissing, "MID"] = newMothers
    }

    # identify connected components and add FAMID column. Keep order unchanged!
    ped = cbind(FAMID = 0, ped)
    connectedIDs = connectedComponentsfb(ped$ID, ped$FID, ped$MID)
    for (i in 1:length(connectedIDs)) ped$FAMID[match(connectedIDs[[i]], ped$ID)] = i


    ### Part2: datamatrix

    if (!is.null(datamatrix)) {
        # sort datamatrix according to ped order
        if (!setequal(familiasped$id, rownames(datamatrix)))
            stop("The datamatrix must have rownames matching the individual ID labels.")
        datamatrix = datamatrix[match(familiasped$id, rownames(datamatrix)), ]

        # convert from factor to character
        allelematrix = sapply(datamatrix, as.character)

        # replace NA with 0
        allelematrix[is.na(allelematrix)] = "0"

        # add empty rows corresponding to new fathers
        addedParents = matrix("0", nrow = nFath + nMoth, ncol = ncol(allelematrix))
        allelematrix = rbind(allelematrix, addedParents)

        # merge with ped columns
        ped = cbind(ped, allelematrix)
    }

    ### Part 3: marker annotations

    annotations = readFamiliasLoci(loci)

    ### Create linkdat object and add plot labels

    x = linkdat(ped, annotations = annotations, verbose = F)
    if (is.linkdat(x))
        x$plot.labels = plotLabels[x$orig.ids] else x = lapply(x, function(xx) {
        xx$plot.labels = plotLabels[xx$orig.ids]
        xx
    })
    x
}

#' @export
#' @rdname Familias2linkdat
readFamiliasLoci = function(loci) {
    if (is.null(loci))
        return(NULL)
    if (inherits(loci, "FamiliasLocus"))
        loci = list(loci)
    annotations = lapply(loci, function(a) {
        malemut = a$maleMutationMatrix
        femalemut = a$femaleMutationMatrix

        if (all(diag(malemut) == 1)) malemut = NULL
        if (all(diag(femalemut) == 1)) femalemut = NULL
        if (is.null(malemut) && is.null(femalemut))
            mutmat = NULL else mutmat = list(male = malemut, female = femalemut)
        list(name = a$locusname, alleles = names(a$alleles), afreq = as.numeric(a$alleles), mutmat = mutmat)
    })
}

#' @export
#' @rdname Familias2linkdat
connectedComponentsfb = function(ID, FID, MID) {
    # Placeholder for final components
    comp = list()

    # Starting point: List of all founders and trios
    temp = lapply(seq_along(ID), function(i) setdiff(c(ID[i], FID[i], MID[i]), 0))  # trios

    while (length(temp) > 0) {
        # Check if first vector overlaps with any of the others
        a = temp[[1]]
        remov = numeric()
        for (j in seq_along(temp)[-1]) if (any(match(a, temp[[j]], nomatch = 0) > 0)) {
            a = unique.default(c(a, temp[[j]]))
            remov = c(remov, j)
        }

        if (length(remov) > 0) {
            # Remove any overlapping vectors, and insert the union as first element
            temp[remov] = NULL
            temp[[1]] = a
        } else {
            # If no overlaps were found, we have a maximal component. Move to comp and remove from
            # temp.
            comp = c(comp, list(sort.default(a)))
            temp[[1]] = NULL
        }
    }
    comp
}

