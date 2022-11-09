#' Modify the pedigree of 'linkdat' objects
#'
#' Functions to modify the pedigree of a 'linkdat' object.
#'
#' When removing an individual, all descendantsfb are also removed as well as
#' founders remaining without offspringfb.
#'
#'
#' @param x A \code{linkdat} object
#' @param father,mother Integers indicating the IDs of parentsfb. If missing, a
#'   new founder individual is created (whose ID will be 1+the largest ID
#'   already in the pedigree).
#' @param noffs A single integer indicating the number of offspringfb to be
#'   created.
#' @param sex,aff Integer vectors indicating the gender and affection statuses
#'   of the offspringfb to be created (recycled if less than \code{noffs}
#'   elements).
#' @param verbose A logical: Verbose output or not.
#' @param keep A character, either 'available' (trimming the pedigree for
#'   unavailable members) or 'affected' (trimming for unaffected members).
#' @param return.ids A logical. If FALSE, the trimmed pedigree is returned as a
#'   new \code{linkdat} object. If TRUE, a vector containing the IDs of
#'   'removable' individuals is returned
#' @param ids individuals
#' @param new a numeric containing new labels to replace those in \code{old}.
#' @param old a numeric containing ID labels to be replaced by those in
#'   \code{new}. If missing, \code{old} is set to \code{x$orig.ids}, i.e. all
#'   members in their original order.
#' @return The modified \code{linkdat} object.
#' @seealso \code{linkdat}, \code{nuclearPed}
#' @import assertthat
#'
#' @name pedModify
NULL


#' @rdname pedModify
#' @export
addOffspring = function(x, father, mother, noffs, ids = NULL, sex = 1, aff = 1, verbose = TRUE) {
    p = as.matrix(x)
    attrs = attributes(p)
    nm = x$nMark
    taken <- oldids <- p[, "ID"]
    if (!missing(father))
        taken = c(taken, father)
    if (!missing(mother))
        taken = c(taken, mother)
    if (!is.null(ids))
        taken = c(taken, ids)
    max_id = max(taken)

    if (missing(father) && missing(mother))
        stop("At least one parent must be an existing pedigree member.")
    if (missing(father))
        father <- max_id <- max_id + 1
    if (missing(mother))
        mother <- max_id <- max_id + 1
    if (any(!is.numeric(father), length(father) != 1))
        stop("Argument 'father' must be a single integer.")
    if (any(!is.numeric(mother), length(mother) != 1))
        stop("Argument 'mother' must be a single integer.")
    if (!any(c(father, mother) %in% oldids))
        stop("At least one parent must be an existing pedigree member.")

    if (missing(noffs) && is.null(ids))
        stop("Number of offspringfb not indicated.")
    if (missing(noffs))
        noffs = length(ids)
    if (is.null(ids))
        ids = (max_id + 1):(max_id + noffs)
    if (length(ids) != noffs)
        stop("Length of 'id' vector must equal number of offspringfb.")
    if (any(ids %in% oldids))
        stop(paste("Individual(s)", ids[ids %in% oldids], "already exist(s)."))

    if (!father %in% oldids) {
        if (verbose)
            cat("Father: Creating new individual with ID", father, "\n")
        p = rbind(p, c(x$famid, father, 0, 0, 1, 1, rep.int(0, nm * 2)))
    }
    if (!mother %in% oldids) {
        if (verbose)
            cat("Mother: Creating new individual with ID", mother, "\n")
        p = rbind(p, c(x$famid, mother, 0, 0, 2, 1, rep.int(0, nm * 2)))
    }
    p = rbind(p, cbind(x$famid, ids, father, mother, sex, aff, matrix(0, ncol = nm * 2, nrow = length(ids))))

    restore_linkdat(p, attrs = attrs)
}


#' @rdname pedModify
#' @export
removeIndividualsfb = function(x, ids, verbose = TRUE) {
    # Remove (one by one) individuals 'ids' and all their descendantsfb. Spouse-founders are
    # removed as well.
    if (any(!ids %in% x$orig.ids))
        stop(paste("Non-existing individuals:", .prettycat(ids[!ids %in% x$orig.ids], "and")))
    pedm = as.matrix(x)

    # Founders without children after 'id' and 'desc' indivs are removed. The redundancy here
    # does not matter.
    desc = numeric(0)
    for (id in ids) {
        desc = c(desc, dd <- descendantsfb(x, id))
        if (verbose)
            cat("Removing", id, if (length(dd) > 0)
                paste("and descendant(s):", .prettycat(dd, "and")), "\n")
    }

    leftover.spousesfb = setdiff(x$orig.ids[x$founders], c(ids, as.numeric(pedm[!x$orig.ids %in%
        c(ids, desc), c("FID", "MID")])))  #founders that are not parentsfb of remaining indivs
    if (verbose && length(leftover.spousesfb) > 0)
        cat("Removing leftover spouse(s):", .prettycat(leftover.spousesfb, "and"), "\n")

    remov = unique(c(ids, desc, leftover.spousesfb))
    restore_linkdat(pedm[-.internalID(x, remov), , drop = F], attrs = attributes(pedm))
}

#' @rdname pedModify
#' @export
trim = function(x, keep = c("available", "affected"), return.ids = FALSE, verbose = TRUE) {
    keep = match.arg(keep)
    if (verbose)
        cat("Trimming pedigree, keeping", keep, "individuals.")
    if (is.singletonfbfb(x) | (keep == "available" && length(x$available) == length(x$orig.ids))) {
        if (verbose)
            cat(" Removed: None\n")
        return(x)
    }
    mysetdiff = function(x, y) x[match(x, y, 0L) == 0L]

    y = linkdat(relabelfb(x$pedigree, x$orig.ids), verbose = F)  # make a copy of x$ped, with original IDs
    y$available = x$available
    while (TRUE) {
        p = y$pedigree
        leavesfb = mysetdiff(p[, "ID"], p[, c("FID", "MID")])
        throw = switch(keep, available = mysetdiff(leavesfb, .internalID(y, y$available)), affected = leavesfb[p[leavesfb,
            "AFF"] != 2])
        if (length(throw) == 0)
            break
        y = removeIndividualsfb(y, y$orig.ids[throw], verbose = FALSE)
    }

    remov = setdiff(x$orig.ids, y$orig.ids)
    if (return.ids)
        return(remov)

    store = as.matrix(x)
    trimmed = store[!store[, "ID"] %in% remov, ]
    if (verbose)
        cat(" Removed:", if (length(remov) > 0)
            .prettycat(remov, "and") else "None", "\n")

    restore_linkdat(trimmed, attrs = attributes(store))
}

#' @rdname pedModify
#' @export
relabelfb = function(x, new, old) {
    islinkdat = is.linkdat(x)
    if (islinkdat) {
        if (length(new) == x$nInd && all(new == x$orig.ids))
            return(x)
        ped = as.matrix(x)
        avail = attr(ped, "available")
    } else ped = x

    orig.ids = ped[, "ID"]
    if (missing(old))
        old = orig.ids
    stopifnot(is.numeric(old), is.numeric(new), length(old) == length(new), !0 %in% new, all(old %in%
        ped[, "ID"]))
    ped[match(old, orig.ids), "ID"] = new

    parentsfb = ped[, c("FID", "MID")]
    ped[, c("FID", "MID")][parentsfb %in% old] <- new[match(parentsfb, old, nomatch = 0)]  #relabelfbing parentsfb

    if (islinkdat) {
        oldavail = avail[avail %in% old]
        avail[avail %in% old] = new[match(oldavail, old)]
        attr(ped, "available") = avail
        return(restore_linkdat(ped))
    } else return(ped)
}


.check_parentsfb_before_children = function(x) {
    if (is.linkdat(x)) {
        ped = x$pedigree # uses internal ordering, but thats OK just for checking
    } else {
        ped = x
        stopifnot(all(c("ID", "FID", "MID") %in% colnames(ped)))
    }
    father_before_child = ped[,'FID'] < ped[,'ID']
    mother_before_child = ped[,'MID'] < ped[,'ID']
    all(father_before_child & mother_before_child)
}


.reorder_parentsfb_before_children = function(x) {
    if(.check_parentsfb_before_children(x))
        return(x)
    if (is.linkdat(x)) {
        ped = as.matrix(x)
        attrs = attributes(ped)
    } else if (all(c("ID", "FID", "MID") %in% colnames(x)))
        ped = x

    N = nrow(ped)
    i = 1
    while (i < N) {
        maxpar = max(match(ped[i, c("FID", "MID")], ped[, "ID"], nomatch = 0))
        if (maxpar > i) {
            ped = ped[c(seq_len(i - 1), (i + 1):maxpar, i, seq_len(N - maxpar) + maxpar), ]
        } else i = i + 1
    }
    if (is.linkdat(x))
        return(restore_linkdat(ped, attrs = attrs))
    else
        return(ped)
}

# Not used?
.merge.linkdat = function(x) {
    # list of linkdats
    if (!is.linkdat.list(x))
        stop("Input must be a list of linkdat objects")
    if (length(x) == 1)
        return(x)
    mnames = lapply(x, function(xx) unlist(lapply(xx$markerdata, attr, "name")))
    common = Reduce(intersect, mnames)
    lapply(x, function(xx) SetMarkersfb(xx, xx$markerdata[getMarkersfb(xx, common)]))
}


#' Functions for modifying availability vectors
#'
#' Functions to set and modify the availability vector of a 'linkdat' object.
#' This vector is used in 'linkage.power' and 'linkageSim', indicating for whom
#' genotypes should be simulated.
#'
#'
#' @param x a \code{linkdat} object
#' @param ids individuals
#' @param available a numeric containing the IDs of available individuals.
#' @return The modified \code{linkdat} object.
#' @seealso \code{plot.linkdat}, \code{linkage.power},
#'   \code{linkageSim}
#'
#'
#' @export
setAvailable = function(x, available) {
    x$available = sort(as.numeric(available))
    x
}

#' @rdname setAvailable
#' @export
swapAvailable = function(x, ids) {
    ava = x$available
    new_ava = c(ava[!ava %in% ids], ids[!ids %in% ava])
    setAvailable(x, new_ava)
}

