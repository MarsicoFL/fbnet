#' Pedigree loops
#'
#' Functions for identifying, breaking and restoring loops in pedigrees.
#'
#' Most of paramlink's handling of pedigree loops is done under the hood - using
#' the functions described here - without need for explicit action from end
#' users. When a linkdat object \code{x} is created, an internal routine detects
#' if the pedigree contains loops, in which case \code{x$hasLoops} is set to
#' TRUE. In analyses of \code{x} where loops must be broken (e.g. lod score
#' computation or marker simulation), this is done automatically by calling
#' \code{breakLoopsfb}.
#'
#' In some cases with complex inbreeding, it can be instructive to plot the
#' pedigree after breaking the loops. Duplicated individuals are plotted with
#' appropriate labels (see examples).
#'
#' The function \code{findLoopBreakersfb} identifies a set of individuals breaking
#' all inbreeding loops, but not marriage loops. These require more machinery
#' for efficient detection, and paramlink does this is a separate function,
#' \code{findLoopBreakersfb2}, utilizing methods from the \code{igraph} package.
#' Since this is rarely needed for most users, \code{igraph} is not imported
#' when loading paramlink, only when \code{findLoopBreakersfb2} is called.
#'
#' In practice, \code{breakLoopsfb} first calls \code{findLoopBreakersfb} and breaks
#' at the returned individuals. If the resulting linkdat object still has loops,
#' \code{findLoopBreakersfb2} is called to break any marriage loops.
#'
#' @param x a \code{\link{linkdat}} object.
#' @param loop_breakers either NULL (resulting in automatic selection of loop
#'   breakers) or a numeric containing IDs of individuals to be used as loop
#'   breakers.
#' @param verbose a logical: Verbose output or not?
#' @return For \code{breakLoopsfb}, a \code{linkdat} object in which the indicated
#'   loop breakers are duplicated. The returned object will also have a non-null
#'   \code{loop_breakers} entry, namely a matrix with the IDs of the original
#'   loop breakers in the first column and the duplicates in the second.
#'
#'   For \code{tieLoopsfb}, a \code{linkdat} object in which any duplicated
#'   individuals (as given in the \code{x$loop_breakers} entry) are merged. For
#'   any linkdat object \code{x}, the call \code{tieLoopsfb(breakLoopsfb(x))} should
#'   return \code{x}.
#'
#'   For \code{pedigreeLoops}, a list containing all inbreeding loops (not
#'   marriage loops) found in the pedigree. Each loop is represented as a list
#'   with elements 'top', a 'bottom' individual, 'pathA' (individuals forming a
#'   path from top to bottom) and 'pathB' (creating a different path from top to
#'   bottom, with no individuals in common with pathA). Note that the number of
#'   loops reported here counts all closed paths in the pedigree and will in
#'   general be larger than the genus of the underlying graph.
#'
#'   For \code{findLoopBreakersfb} and \code{findLoopBreakersfb2}, a numeric vector
#'   of individual ID's.
#'
#'
#' @export
pedigreeLoops = function(x) {
    dls = .descentPaths(x, 1:x$nInd, original.ids = FALSE)

    loops = list()
    for (id in 1:x$nInd) {
        if (length(dl <- dls[[id]]) == 1)
            next
        pairs = .comb2(length(dl))
        for (p in 1:nrow(pairs)) {
            p1 = dl[[pairs[p, 1]]]
            p2 = dl[[pairs[p, 2]]]
            if (p1[2] == p2[2])
                next
            inters = p1[match(p2, p1, 0L)][-1]  #intersecting path members, excluding id
            if (length(inters) == 0)
                next else {
                top = x$orig.ids[p1[1]]
                bottom = x$orig.ids[inters[1]]
                pathA = p1[seq_len(which(p1 == inters[1]) - 2) + 1]  #without top and bottom. Seq_len to deal with the 1:0 problem.
                pathB = p2[seq_len(which(p2 == inters[1]) - 2) + 1]
                loops = c(loops, list(list(top = top, bottom = bottom, pathA = x$orig.ids[pathA],
                  pathB = x$orig.ids[pathB])))
            }
        }
    }
    unique(loops)
}

#' @export
#' @rdname pedigreeLoops
breakLoopsfb = function(x, loop_breakers = NULL, verbose = TRUE) {
    if (is.singletonfbfb(x))
        stop("This function does not apply to singleton objects.")
    automatic = is.null(loop_breakers)
    if (automatic) {
        if (!x$hasLoops)
            return(x)
        loop_breakers = findLoopBreakersfb(x)
        if (length(loop_breakers) == 0) {
            if (verbose)
                cat("Marriage loops detected, trying different selection method.\n")
            loop_breakers = findLoopBreakersfb2(x)
        }
    }

    if (length(loop_breakers) == 0)
        stop("Loop breaking unsuccessful.")

    if (any(loop_breakers %in% x$orig.ids[x$founders]))
        stop("Pedigree founders cannot be loop breakers.")

    if (verbose)
        cat(ifelse(automatic, "Selected", "User specified"), "loop breakers:", loop_breakers,
            "\n")

    pedm = as.matrix(x)  #data.frame(x, famid=T, missing=0)
    attrs = attributes(pedm)  #all attributes except 'dim'
    dup_pairs = x$loop_breakers  #normally = NULL at this stage
    for (id in loop_breakers) {
        dup_id = max(pedm[, "ID"]) + 1
        dup_pairs = rbind(dup_pairs, c(id, dup_id))
        intern = match(id, pedm[, "ID"])  #don't use .internalID here, since pedm changes all the time
        sex_col = pedm[intern, "SEX"] + 2  # FID column if 'intern' is male; MID if female

        pedm = pedm[c(1:intern, intern, (intern + 1):nrow(pedm)), ]
        pedm[intern + 1, c("ID", "FID", "MID")] = c(dup_id, 0, 0)
        pedm[pedm[, sex_col] == id, sex_col] = dup_id
    }
    newx = restore_linkdat(pedm, attrs = attrs)
    newx$loop_breakers = dup_pairs
    if (automatic)
        return(breakLoopsfb(newx, verbose = verbose))
    newx
}

#' @export
#' @rdname pedigreeLoops
tieLoopsfb = function(x) {
    dups = x$loop_breakers
    if (is.null(dups) || nrow(dups) == 0) {
        cat("No loops to tie\n")
        return(x)
    }
    if (!all(dups %in% x$orig.ids))
        stop("Something's wrong: Duplicated individuals no longer in pedigree.")
    pedm = as.matrix(x)
    attrs = attributes(pedm)

    origs = dups[, 1]
    copies = dups[, 2]
    pedm = pedm[-match(copies, pedm[, "ID"]), ]
    for (i in 1:length(origs)) {
        orig = origs[i]
        copy = copies[i]
        sex = pedm[pedm[, "ID"] == orig, "SEX"]
        pedm[pedm[, sex + 2] == copy, sex + 2] = orig
    }
    restore_linkdat(pedm, attrs = attrs)
}

#' @export
#' @rdname pedigreeLoops
findLoopBreakersfb = function(x) {
    loopdata = pedigreeLoops(x)
    # write each loop as vector exluding top/bottom
    loops = lapply(loopdata, function(lo) c(lo$pathA, lo$pathB))
    bestbreakers = numeric()
    while (length(loops) > 0) {
        # add the individual occuring in most loops
        best = which.max(tabulate(unlist(loops)))
        bestbreakers = c(bestbreakers, best)
        loops = loops[sapply(loops, function(vec) !best %in% vec)]
    }
    bestbreakers
}

#' @export
#' @rdname pedigreeLoops
findLoopBreakersfb2 = function(x) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        cat("This pedigree has marriage loops. The package 'igraph' must be installed for automatic selection of loop breakers.")
        return(numeric())
    }

    ids = x$orig.ids
    N = max(ids)
    nonf = ids[x$nonfounders]
    breakers = numeric()

    ped2edge = function(p) {
        # input: ped-matrise UTEN founder-rader
        pp = cbind(p, paste(p[, "FID"], p[, "MID"], sep = "+"))
        edge.children = pp[, c(4, 1), drop = F]
        edge.marriage = unique.matrix(rbind(pp[, c(2, 4), drop = F], pp[, c(3, 4), drop = F]))
        rbind(edge.marriage, edge.children)
    }

    p = as.matrix(x, keep = F)[x$nonfounders, c("ID", "FID", "MID"), drop = F]
    while (TRUE) {
        g = igraph::graph_from_edgelist(ped2edge(p))
        loop = igraph::girth(g)$circle
        if (length(loop) == 0)
            break
        good.breakers = intersect(loop$name, nonf)
        if (length(good.breakers) == 0)
            stop("This pedigree requires founders as loop breakers, which is not implemented in paramlink yet. Sorry!")
        b = good.breakers[1]
        breakers = c(breakers, b)

        b.is.parent = p == b & col(p) > 1
        p[b.is.parent] = as.character(N <- N + 1)
    }
    breakers
}

