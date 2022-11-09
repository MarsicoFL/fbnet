#' Pedigree likelihood
#'
#' Calculates various forms of pedigree likelihoods based on paramlink functions.
#'
#' All likelihoods are calculated using the Elston-Stewart algorithm.
#'
#'
#' @param x a \code{linkdat} object, a singleton object, or
#'   a list of such objects. In \code{likelihood_LINKAGE}, \code{x} must be a
#'   \code{linkdat} object, with \code{x$model} different from NULL.
#' @param locus1 a \code{marker} object compatible with \code{x}. If
#'   \code{x} is a list, then \code{locus1} must be a list of corresponding
#'   \code{marker} objects.
#' @param locus2 either NULL, the character 'disease', or a markerfb
#'   object compatible with \code{x}. See Details.
#' @param theta the recombination rate between locus1 and locus2 (in
#'   \code{likelihood_LINKAGE}: between the marker and the disease locus).  To
#'   make biological sense theta should be between 0 and 0.5.
#' @param eliminate mostly for internal use: a non-negative integer indicating
#'   the number of iterations in the internal genotype-compatibility algorithm.
#'   Positive values can save time if \code{partialmarker} is non-empty and the
#'   number of alleles is large.
#' @param logbase a numeric, or NULL. If numeric the log-likelihood is returned,
#'   with \code{logbase} as basis for the logarithm.
#' @param loop_breakers a numeric containing IDs of individuals to be used as
#'   loop breakers. If NULL, automatic selection of loop breakers will be
#'   performed. 
#' @param startdata for internal use.
#'   linkage computations with few-allelic markers.
#' @param \dots further arguments.
#' @return The likelihood of the data. If the parameter \code{logbase} is a
#'   positive number, the output is \code{log(likelihood, logbase)}.
#'

#' @export
#' @rdname likelihood
likelihood.linkdat = function(x, locus1, locus2 = NULL, theta = NULL, startdata = NULL, eliminate = 0,
    logbase = NULL, loop_breakers = NULL, ...) {
    if (!inherits(locus1, "marker") && is.numeric(locus1))
        locus1 = x$markerdata[[locus1]]
    if (!inherits(locus2, "marker") && is.numeric(locus2))
        locus2 = x$markerdata[[locus2]]

    # analysisType = switch(class(locus2), `NULL` = 1, character = 2, marker = 3)
    if(inherits(locus2, "marker"))
      analysisType = 3
    else if(is.character(locus2)) {
      analysisType = 2
      locus2 = NULL  # 'disease'
    }
    else
      analysisType = 1

    locus1 = .reduce_alleles(locus1)
    locus2 = .reduce_alleles(locus2)  # unchanged if NULL

    if (x$hasLoops) {
        cat("Tip: To optimize speed consider breaking loops before calling 'likelihood'. See ?breakLoopsfb.\n")
        m = list(locus1)
        m[[2]] = locus2  # no effect if NULL
        x = breakLoopsfb(SetMarkersfb(x, m), loop_breakers = loop_breakers, verbose = TRUE)
        locus1 = x$markerdata[[1]]
        if (analysisType == 3)
            locus2 = x$markerdata[[2]]
    }

    chrom = if (identical(attr(locus1, "chrom"), 23))
        "X" else "AUTOSOMAL"
    SEX = x$pedigree[, "SEX"]
    mutmat = attr(locus1, "mutmat")
    inform_subnucs = x$subnucs

    if (is.null(startdata)) {
        if (analysisType != 2) {
            inform = .informative(x, locus1, locus2)
            inform_subnucs = inform$subnucs
            x$founders = c(x$founders, inform$newfounders)
            x$nonfounders = .mysetdiff(x$nonfounders, inform$newfounders)
        }
        dat = switch(analysisType, .startdata_M(x, marker = locus1, eliminate = eliminate),
            .startdata_MD(x, marker = locus1, eliminate = eliminate), .startdata_MM(x, marker1 = locus1,
                marker2 = locus2, eliminate = eliminate))
    } else dat = startdata

    peelFUN = switch(analysisType, function(dat, sub) .peel_M(dat, sub, chrom, SEX, mutmat = mutmat),
        function(dat, sub) .peel_MD(dat, sub, theta, chrom, SEX), function(dat, sub) .peel_MM(dat,
            sub, theta, chrom, SEX))

    if (attr(dat, "impossible"))
        return(ifelse(is.numeric(logbase), -Inf, 0))

    if (is.null(dups <- x$loop_breakers)) {
        for (sub in inform_subnucs) {
            dat = peelFUN(dat, sub)
            if (sub$pivtype > 0 && attr(dat, "impossible"))
                return(ifelse(is.numeric(logbase), -Inf, 0))
        }
        likelihood = dat
    } else {
        two2one = function(matr) if (is.matrix(matr))
            1000 * matr[1, ] + matr[2, ] else matr  #if input is vector (i.e. X-linked male genotypes), return it unchanged
        origs = match(dups[, 1], x$orig.ids)
        copies = match(dups[, 2], x$orig.ids)

        # For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy.
        # Then take cross product of these vectors.
        loopgrid = fast.grid(lapply(seq_along(origs), function(i) {
            ori = two2one(dat[[c(origs[i], 1)]])
            seq_along(ori)[ori %in% two2one(dat[[c(copies[i], 1)]])]
        }), as.list = TRUE)

        likelihood = 0
        for (r in loopgrid) {
            # r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
            dat1 = dat
            attr(dat1, "impossible") = FALSE
            for (i in seq_along(origs)) {
                orig.int = origs[i]
                copy.int = copies[i]
                orighap = dat[[orig.int]]$hap
                origprob = dat[[orig.int]]$prob
                hap = if (is.matrix(orighap))
                  orighap[, r[i], drop = F] else orighap[r[i]]
                prob = origprob[r[i]]
                if (sum(prob) == 0)
                  print("Loop-loekke: Alle sannsynligheter er null. Magnus lurer paa om dette kan gi feilmelding.")
                dat1[[orig.int]] = list(hap = hap, prob = prob)
                dat1[[copy.int]] = list(hap = hap, prob = 1)
            }
            for (sub in inform_subnucs) {
                pivtype = sub$pivtype
                dat1 = peelFUN(dat1, sub)
                if (pivtype > 0 && attr(dat1, "impossible"))
                  {
                    break
                  }  #if impossible data - break out of ES-algorithm and go to next r in loopgrid.
                if (pivtype == 0)
                  likelihood = likelihood + dat1
            }
        }
    }
    if (is.numeric(logbase))
        log(likelihood, logbase) else likelihood
}

#' @export
#' @rdname likelihood
likelihood.singleton = function(x, locus1, logbase = NULL, ...) {
    if (!inherits(locus1, "marker") && is.numeric(locus1))
        locus1 = x$markerdata[[locus1]]
    if (is.null(locus1) || all(locus1 == 0))
        return(if (is.numeric(logbase)) 0 else 1)

    m = locus1
    chrom = as.integer(attr(m, "chrom"))
    afreq = attr(m, "afreq")
    if (identical(chrom, 23L) && x$pedigree[, "SEX"] == 1) {
        # X chrom and male
        if (all(m > 0) && m[1] != m[2])
            stop("Heterozygous genotype detected for X-linked marker in male individual.")
        res = afreq[m[1]]
    } else if (0 %in% m) {
        p = afreq[m[m != 0]]
        res = p^2 + 2 * p * (1 - p)
    } else res = prod(afreq[m]) * ifelse(m[1] != m[2], 2, 1)  # assumes HWE
    return(if (is.numeric(logbase)) log(res, logbase) else res)
}

# #' @export
# #' @rdname likelihood
# likelihood.list = function(x, locus1, locus2 = NULL, ..., returnprod = TRUE) {
#     if (!is.linkdat.list(x))
#        stop("x must be either a 'linkdat' object, a 'singleton' object, or a list of such")
#    if (is.atomic(locus1))
#        locus1 = rep(list(locus1), length = length(x))
#    if (is.atomic(locus2))
#        locus2 = rep(list(locus2), length = length(x))  # Note: NULL is atomic
#    liks = vapply(1:length(x), function(i) likelihood(x[[i]], locus1[[i]], locus2[[i]], ...),
#        numeric(1))
#
#    if (returnprod)
#        return(prod(liks)) else liks
#}


#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

.startdata_M = function(x, marker, eliminate = 0) {
    afreq = attr(marker, "afreq")
    chromX = identical(attr(marker, "chrom"), 23)

    impossible = FALSE
    if (chromX) {
        glist = .build_genolist_X(x, marker, eliminate)
        if (attr(glist, "impossible")) {
            dat = list()
            attr(dat, "impossible") = TRUE
            return(dat)
        }
        sex = x$pedigree[, "SEX"]
        dat = lapply(1:x$nInd, function(i) {
            h = glist[[i]]
            if (i %in% x$founders) {
                prob = switch(sex[i], afreq[h], afreq[h[1, ]] * afreq[h[2, ]] * ((h[1, ] !=
                  h[2, ]) + 1))
                if (sum(prob) == 0)
                  impossible = TRUE
            } else prob = rep.int(1, length(h)/sex[i])
            list(hap = h, prob = as.numeric(prob))
        })
    } else {
        glist = .build_genolist(x, marker, eliminate)
        if (attr(glist, "impossible")) {
            dat = list()
            attr(dat, "impossible") = TRUE
            return(dat)
        }
        dat = lapply(1:x$nInd, function(i) {
            h = glist[[i]]
            if (i %in% x$founders) {
                prob = afreq[h[1, ]] * afreq[h[2, ]] * ((h[1, ] != h[2, ]) + 1)
                if (sum(prob) == 0)
                  impossible = TRUE
            } else prob = rep.int(1, ncol(h))
            list(hap = h, prob = as.numeric(prob))
        })
    }
    attr(dat, "impossible") = impossible
    dat
}


.startdata_MD = function(x, marker, eliminate = 0) {
    startprob <- function(h, model, afreq, aff, founder) {
        pat = h[1, ]
        mat = h[2, ]
        al1 = abs(pat)
        al2 = abs(mat)
        d.no = (pat < 0) + (mat < 0)
        prob = switch(aff + 1, rep.int(1, length(d.no)), (1 - model$penetrances)[d.no + 1],
            model$penetrances[d.no + 1])
        if (founder)
            prob = prob * afreq[al1] * afreq[al2] * ((al1 != al2) + 1) * model$dfreq^d.no *
                (1 - model$dfreq)^(2 - d.no)
        as.numeric(prob)
    }

    startprob_X <- function(h, model, afreq, sex, aff, founder) {
        switch(sex, {
            mat = h  #vector
            d.no = as.numeric(mat < 0)
            prob <- switch(aff + 1, rep.int(1, length(d.no)), (1 - model$penetrances$male)[d.no +
                1], model$penetrances$male[d.no + 1])
            if (founder) prob <- prob * afreq[abs(mat)] * c(1 - model$dfreq, model$dfreq)[d.no +
                1]
        }, {
            pat = h[1, ]
            mat = h[2, ]
            al1 = abs(pat)
            al2 = abs(mat)
            d.no = (pat < 0) + (mat < 0)
            prob <- switch(aff + 1, rep.int(1, length(d.no)), (1 - model$penetrances$female)[d.no +
                1], model$penetrances$female[d.no + 1])
            if (founder) prob <- prob * afreq[al1] * afreq[al2] * ((al1 != al2) + 1) * model$dfreq^d.no *
                (1 - model$dfreq)^(2 - d.no)
        })
        as.numeric(prob)
    }

    ## MAIN ###
    afreq = attr(marker, "afreq")
    chromX = identical(attr(marker, "chrom"), 23)
    AFF = x$pedigree[, "AFF"]
    FOU = (1:x$nInd) %in% x$founders
    dlist = cbind(c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))  #D=-1; N=1
    impossible = FALSE
    if (chromX) {
        glist = .build_genolist_X(x, marker, eliminate)
        if (attr(glist, "impossible")) {
            dat = list()
            attr(dat, "impossible") = TRUE
            return(dat)
        }

        SEX = x$pedigree[, "SEX"]
        dat = lapply(1:x$nInd, function(i) {
            switch(SEX[i], {
                hap = c(glist[[i]], -glist[[i]])
                prob = startprob_X(hap, model = x$model, afreq = afreq, sex = 1, aff = AFF[i],
                  founder = FOU[i])
                list(hap = hap[prob > 0], prob = prob[prob > 0])
            }, {
                gl = ncol(glist[[i]])
                hap = glist[[i]][, rep(1:gl, each = 4), drop = F] * dlist[, rep(1:4, times = gl),
                  drop = F]
                prob = startprob_X(hap, model = x$model, afreq = afreq, sex = 2, aff = AFF[i],
                  founder = FOU[i])
                if (sum(prob) == 0) impossible = TRUE
                list(hap = hap[, prob > 0, drop = F], prob = as.numeric(prob[prob > 0]))
            })
        })
    } else {
        glist = .build_genolist(x, marker, eliminate)
        if (attr(glist, "impossible")) {
            dat = list()
            attr(dat, "impossible") = TRUE
            return(dat)
        }

        dat = lapply(1:x$nInd, function(i) {
            gl = ncol(glist[[i]])
            hap = glist[[i]][, rep(1:gl, each = 4), drop = F] * dlist[, rep(1:4, times = gl),
                drop = F]
            prob = startprob(hap, model = x$model, afreq = afreq, aff = AFF[i], founder = FOU[i])
            if (sum(prob) == 0)
                impossible = TRUE
            list(hap = hap[, prob > 0, drop = F], prob = as.numeric(prob[prob > 0]))
        })
    }
    attr(dat, "impossible") = impossible
    dat
}


.startdata_MM = function(x, marker1, marker2, eliminate = 0) {

    startprob <- function(h, afreq1, afreq2, founder) {
        if (founder) {
            m1_1 = h[1, ]
            m1_2 = h[2, ]
            m2_1 = h[3, ]
            m2_2 = h[4, ]
            hetfact = ((m1_1 != m1_2 | m2_1 != m2_2) + 1)  #multiply with two if heteroz for at least 1 marker. If heteroz for both, then both phases are included in h, hence the factor 2 (not 4) in this case as well.
            prob = afreq1[m1_1] * afreq1[m1_2] * afreq2[m2_1] * afreq2[m2_2] * hetfact
            return(as.numeric(prob))
        }
        return(rep.int(1, ncol(h)))
    }

    startprob_X <- function(h, afreq1, afreq2, sex, founder) {
        if (founder && sex == 1)
            return(as.numeric(afreq1[h[1, ]] * afreq2[h[2, ]]))
        startprob(h, afreq1, afreq2, founder)
    }

    afreq1 = attr(marker1, "afreq")
    afreq2 = attr(marker2, "afreq")
    chromX = identical(attr(marker1, "chrom"), 23)
    impossible = FALSE

    if (chromX) {
        sex = x$pedigree[, "SEX"]
        m1_list = .build_genolist_X(x, marker1, eliminate)
        m2_list = .build_genolist_X(x, marker2, eliminate)
        if (attr(m1_list, "impossible") || attr(m2_list, "impossible")) {
            dat = list()
            attr(dat, "impossible") = TRUE
            return(dat)
        }

        dat = lapply(1:x$nInd, function(i) {
            sexi = sex[i]
            founder = i %in% x$founders
            h1 = m1_list[[i]]
            h2 = m2_list[[i]]
            if (sexi == 1)
                hap = rbind(rep(h1, each = length(h2)), rep(h2, times = length(h1)))  #matrix with two rows: m1, m2
 else {
                hl1 = dim(h1)[2]
                hl2 = dim(h2)[2]
                hap = rbind(h1[, rep(seq_len(hl1), each = hl2), drop = F], h2[, rep(seq_len(hl2),
                  , times = hl1), drop = FALSE])  #matrix with four rows: m1_1, m1_2, m2_1, m2_2
                if (founder) {
                  # doubly heterozygous founders: Include the other phase as well. (This is necessary since
                  # .build_genolist returns unordered genotypes for founders.)
                  doublyhet = hap[1, ] != hap[2, ] & hap[3, ] != hap[4, ]
                  if (any(doublyhet))
                    hap = cbind(hap, hap[c(1, 2, 4, 3), doublyhet, drop = FALSE])
                }
            }
            prob = startprob_X(hap, afreq1 = afreq1, afreq2 = afreq2, sex = sexi, founder = founder)
            keep = prob > 0
            if (!any(keep))
                impossible = TRUE
            list(hap = hap[, keep, drop = F], prob = as.numeric(prob[keep]))
        })
    } else {
        m1_list = .build_genolist(x, marker1, eliminate)
        m2_list = .build_genolist(x, marker2, eliminate)
        if (attr(m1_list, "impossible") || attr(m2_list, "impossible")) {
            dat = list()
            attr(dat, "impossible") = TRUE
            return(dat)
        }

        dat = lapply(1:x$nInd, function(i) {
            h1 = m1_list[[i]]
            hl1 = dim(h1)[2]
            h2 = m2_list[[i]]
            hl2 = dim(h2)[2]
            hap = rbind(h1[, rep(seq_len(hl1), each = hl2), drop = F], h2[, rep(seq_len(hl2),
                , times = hl1), drop = FALSE])  #matrix with four rows: m1_1, m1_2, m2_1, m2_2
            if (i %in% x$founders) {
                # doubly heterozygous founders: Include the other phase as well. (This is necessary since
                # .build_genolist returns unordered genotypes for founders.)
                doublyhet = hap[1, ] != hap[2, ] & hap[3, ] != hap[4, ]
                if (any(doublyhet))
                  hap = cbind(hap, hap[c(1, 2, 4, 3), doublyhet, drop = FALSE])
            }
            prob = startprob(hap, afreq1 = afreq1, afreq2 = afreq2, founder = (i %in% x$founders))
            keep = prob > 0
            if (!any(keep))
                impossible = TRUE
            list(hap = hap[, keep, drop = F], prob = as.numeric(prob[keep]))
        })
    }
    attr(dat, "impossible") = impossible
    dat
}


#### .BUILD_GENOLIST and ELIMINATE

.build_genolist <- function(x, marker, eliminate = 0) {
    # mm: marker matrix, dim = (nInd , 2)
    n = attr(marker, "nalleles")
    nseq = seq_len(n)

    .COMPLETE = {
        tmp1 = rep(nseq, each = n)
        tmp2 = rep.int(nseq, times = n)
        fath = c(tmp1, tmp2)
        moth = c(tmp2, tmp1)
        rbind(fath, moth, deparse.level = 0)[, !duplicated.default(fath * 1000 + moth), drop = F]  #faster than unique(m, MARGIN=2)
    }
    genolist = lapply(1:x$nInd, function(i) {
        g_i = marker[i, ]
        m = switch(sum(g_i == 0) + 1, {
            a = g_i[1]
            b = g_i[2]
            if (a == b) cbind(g_i, deparse.level = 0) else cbind(g_i, c(b, a), deparse.level = 0)
        }, {
            nz = g_i[g_i != 0]
            rbind(c(nseq, rep.int(nz, n - 1)), c(rep.int(nz, n), nseq[-nz]), deparse.level = 0)
        }, {
            .COMPLETE
        })
        if ((i %in% x$founders) && (!x$orig.ids[i] %in% x$loop_breakers[, 2]))
            m = m[, m[1, ] <= m[2, ], drop = FALSE]
        m
    })
    attr(genolist, "impossible") = FALSE

    # If mutations, don't eliminate any genotypes
    if (!is.null(attr(marker, "mutmat")))
        return(genolist)

    .eliminate(x, genolist, n, repeats = eliminate)
}


.eliminate = function(x, genolist, nall, repeats = 0) {
    if (repeats == 0 || attr(genolist, "impossible"))
        return(genolist)
    offs = lapply(1:x$nInd, function(i) offspringfb(x, i, original.id = FALSE))
    ncols_ny = unlist(lapply(genolist, ncol))
    p = x$pedigree
    informative = logical(x$nInd)
    for (k in seq_len(repeats)) {
        ncols = ncols_ny
        informative[x$founders] = (ncols[x$founders] < nall * (nall + 1)/2)
        informative[x$nonfounders] = (ncols[x$nonfounders] < nall^2)
        for (i in 1:x$nInd) {
            if (ncols[i] == 1)
                next
            g = genolist[[i]]
            kjonn = p[i, "SEX"]
            if (i %in% x$nonfounders && informative[far <- p[i, "FID"]])
                g = g[, g[1, ] %in% genolist[[far]][1, ] | g[1, ] %in% genolist[[far]][2, ],
                  drop = F]
            if (i %in% x$nonfounders && informative[mor <- p[i, "MID"]])
                g = g[, g[2, ] %in% genolist[[mor]][1, ] | g[2, ] %in% genolist[[mor]][2, ],
                  drop = F]
            barn = offs[[i]]
            for (b in barn[informative[barn]]) {
                g = g[, g[1, ] %in% genolist[[b]][kjonn, ] | g[2, ] %in% genolist[[b]][kjonn,
                  ], drop = F]
            }
            genolist[[i]] = g
        }
        ncols_ny = unlist(lapply(genolist, ncol))
        if (any(ncols_ny == 0)) {
            attr(genolist, "impossible") = TRUE
            return(genolist)
        }
        if (sum(ncols_ny) == sum(ncols))
            return(genolist)
    }
    genolist
}


.reduce_alleles = function(marker) {
    if (all(marker != 0))
        return(marker)  # no reduction needed (OK!)
    attrs = attributes(marker)

    if (!is.null(attrs$mutmat)) {
        malem = attrs$mutmat$male
        femalem = attrs$mutmat$female
        male_lump = identical(attr(malem, "lumpability"), "always")
        female_lump = identical(attr(femalem, "lumpability"), "always")
        if (!male_lump || !female_lump)
            return(marker)
    }
    orig_alleles = attrs$alleles

    # indices of observed alleles
    present = sort.int(setdiff(unique.default(as.numeric(marker)), 0))
    if (length(present) >= length(orig_alleles) - 1)
        return(marker)  # return unchanged if all, or all but one, are observed

    redund = setdiff(1:attrs$nalleles, present)
    dummylab = paste(orig_alleles[redund], collapse = "_")

    if (length(present) == 0) {
        new_marker = rep.int(0, length(marker))
        attributes(new_marker) = modifyList(attrs, list(alleles = dummylab, nalleles = 1, afreq = 1))
        if (!is.null(attrs$mutmat)) {
            mm = matrix(1, dimnames = list(dummylab, dummylab))
            attr(new_marker, "mutmat") = list(male = mm, female = mm)
        }
        return(new_marker)
    }

    new_marker = match(marker, present, nomatch = 0)
    new_alleles = c(orig_alleles[present], dummylab)
    present_freq = attrs$afreq[present]
    new_freq = c(present_freq, 1 - sum(present_freq))
    n = length(present) + 1

    attributes(new_marker) = modifyList(attrs, list(alleles = new_alleles, nalleles = n, afreq = new_freq))

    if (!is.null(attrs$mutmat)) {
        if (male_lump) {
            mm = malem[c(present, redund[1]), c(present, redund[1])]
            mm[, n] = 1 - rowSums(mm[, -n, drop = F])
        }
        if (female_lump) {
            mf = femalem[c(present, redund[1]), c(present, redund[1])]
            mf[, n] = 1 - rowSums(mf[, -n, drop = F])
        }
        # for(i in 1:(n-1)) { present_i = present[i] #m_weight = f_weight = attrs$afreq[redund]
        # m_weight = malem[present_i, redund] f_weight = femalem[present_i, redund] mm[n, i] =
        # (m_weight/sum(m_weight)) %*% malem[redund, present_i] mf[n, i] = (f_weight/sum(f_weight))
        # %*% femalem[redund, present_i] }
        attr(new_marker, "mutmat") = list(male = mm, female = mf)
    }
    new_marker
}


#------------X-linked-------------------

.build_genolist_X <- function(x, marker, eliminate) {
    # marker: marker matrix, dim = (nInd , 2).
    n = attr(marker, "nalleles")
    nseq = seq_len(n)

    .COMPLETE = {
        tmp1 = rep(nseq, each = n)
        tmp2 = rep.int(nseq, times = n)
        fath = c(tmp1, tmp2)
        moth = c(tmp2, tmp1)
        rbind(fath, moth, deparse.level = 0)[, !duplicated.default(fath * 1000 + moth), drop = F]  #faster than unique(m, MARGIN=2)
    }
    SEX = x$pedigree[, "SEX"]
    females = (1:x$nInd)[SEX == 2]
    genolist = list()
    genolist[SEX == 1 & marker[, 1] == 0] = list(nseq)
    genolist[SEX == 1 & marker[, 1] != 0] = marker[SEX == 1 & marker[, 1] != 0, 1]
    genolist[females] = lapply(females, function(i) {
        g_i = marker[i, ]
        m = switch(sum(g_i == 0) + 1, {
            a = g_i[1]
            b = g_i[2]
            if (a == b) cbind(g_i, deparse.level = 0) else cbind(g_i, c(b, a), deparse.level = 0)
        }, {
            nz = g_i[g_i != 0]
            rbind(c(nseq, rep.int(nz, n - 1)), c(rep.int(nz, n), nseq[-nz]), deparse.level = 0)
        }, {
            .COMPLETE
        })
        if ((i %in% x$founders) && (!x$orig.ids[i] %in% x$loop_breakers))
            m = m[, m[1, ] <= m[2, ], drop = FALSE]
        m
    })
    attr(genolist, "impossible") = FALSE

    # If mutations, don't eliminate any genotypes
    if (!is.null(attr(marker, "mutmat")))
        return(genolist)

    .eliminate_X(x, genolist, n, eliminate)
}


.eliminate_X = function(x, genolist, nall, repeats = 0) {
    if (repeats == 0 || attr(genolist, "impossible"))
        return(genolist)
    SEX = x$pedigree[, "SEX"]
    FID = x$pedigree[, "FID"]
    MID = x$pedigree[, "MID"]
    males = (1:x$nInd)[SEX == 1]
    females = (1:x$nInd)[SEX == 2]
    fem_fou = .myintersect(females, x$founders)
    fem_nonfou = .myintersect(females, x$nonfounders)
    offs = lapply(1:x$nInd, function(i) offspringfb(x, i, original.id = FALSE))
    p = x$pedigree
    informative = logical(x$nInd)
    ncols_ny = lengths(genolist)/SEX  #males are vectors, females matrices w/ 2 rows
    for (k in seq_len(repeats)) {
        ncols = ncols_ny
        informative[males] = (ncols[males] < nall)
        informative[fem_fou] = (ncols[fem_fou] < nall * (nall + 1)/2)
        informative[fem_nonfou] = (ncols[fem_nonfou] < nall^2)
        for (i in males) {
            if (ncols[i] == 1)
                next
            g = genolist[[i]]
            if (i %in% x$nonfounders && informative[mor <- MID[i]])
                g = g[g %in% genolist[[mor]][1, ] | g %in% genolist[[mor]][2, ]]
            barn = offs[[i]]
            for (b in barn[informative[barn] & SEX[barn] == 2]) g = g[g %in% genolist[[b]][1,
                ]]
            genolist[[i]] = g
        }
        for (i in females) {
            if (ncols[i] == 1)
                next
            g = genolist[[i]]
            if (i %in% x$nonfounders && informative[far <- FID[i]])
                g = g[, g[1, ] %in% genolist[[far]], drop = F]
            if (i %in% x$nonfounders && informative[mor <- MID[i]])
                g = g[, g[2, ] %in% genolist[[mor]][1, ] | g[2, ] %in% genolist[[mor]][2, ],
                  drop = F]
            barn = offs[[i]]
            for (b in barn[informative[barn]]) {
                if (SEX[b] == 1)
                  g = g[, g[1, ] %in% genolist[[b]] | g[2, ] %in% genolist[[b]], drop = F] else g = g[, g[1, ] %in% genolist[[b]][2, ] | g[2, ] %in% genolist[[b]][2,
                  ], drop = F]
            }
            genolist[[i]] = g
        }
        ncols_ny = lengths(genolist)/SEX
        if (any(ncols_ny == 0)) {
            attr(genolist, "impossible") = TRUE
            return(genolist)
        }
        if (sum(ncols_ny) == sum(ncols))
            return(genolist)
    }
    genolist
}


#### PEELING FUNCTIONS (one for each case: single Marker, Marker-Disease, Marker-Marker)

.peel_M <- function(dat, sub, chrom, SEX, mutmat = NULL) {
    far = sub[["father"]]
    mor = sub[["mother"]]
    offs = sub[["offspringfb"]]
    piv = sub[["pivot"]]
    pivtype = sub[["pivtype"]]
    if (pivtype == 3)
        offs = offs[offs != piv]  #pivtype indicates who is pivot: 0 = none; 1 = father; 2 = mother; 3 = an offspring
    farh = dat[[c(far, 1)]]
    morh = dat[[c(mor, 1)]]
    likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
    dims = dim(likel)
    fa_len = dims[1L]
    mo_len = dims[2L]

    .trans_M <- function(parent, childhap, mutmat = NULL) {
        # parent = matrix with 2 rows; childhap = vector of any length (parental allele); mutmat =
        # mutation matrix
        if (is.null(mutmat))
            unlist(lapply(seq_len(ncol(parent)), function(i) ((parent[1, i] == childhap) +
                (parent[2, i] == childhap))/2)) else unlist(lapply(seq_len(ncol(parent)), function(i) (mutmat[parent[1, i], childhap] +
            mutmat[parent[2, i], childhap])/2))
    }

    switch(chrom, AUTOSOMAL = {
        for (datb in dat[offs]) {
            bh = datb[[1]]
            bp = datb[[2]]
            bl = length(bp)
            trans_pats = .trans_M(farh, bh[1, ], mutmat = mutmat$male)
            trans_mats = .trans_M(morh, bh[2, ], mutmat = mutmat$female)
            dim(trans_mats) = c(bl, mo_len)
            trans_mats_rep = as.numeric(do.call(rbind, rep(list(trans_mats), fa_len)))
            mm = .colSums((trans_pats * bp) * trans_mats_rep, bl, fa_len * mo_len)
            likel = likel * mm
        }
        if (pivtype == 0) return(sum(likel))

        res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
            {
                pivh = dat[[c(piv, 1)]]
                pivp = dat[[c(piv, 2)]]
                pi_len = length(pivp)
                T = numeric(fa_len * mo_len * pi_len)
                dim(T) = c(fa_len, mo_len, pi_len)
                trans_pats = .trans_M(farh, pivh[1, ], mutmat = mutmat$male)
                dim(trans_pats) = c(pi_len, fa_len)
                trans_mats = .trans_M(morh, pivh[2, ], mutmat = mutmat$female)
                dim(trans_mats) = c(pi_len, mo_len)
                for (i in seq_len(fa_len)) {
                  transpat = trans_pats[, i]
                  for (j in seq_len(mo_len)) T[i, j, ] = transpat * trans_mats[, j]
                }
                arr = as.vector(T) * as.vector(likel)
                dim(arr) = dim(T)
                res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[piv]]
                res * pivp
            })
        pivhap_update = dat[[c(piv, 1)]][, res > 0, drop = F]
    }, X = {
        for (b in offs) {
            datb = dat[[b]]
            bh = datb[[1]]
            bp = datb[[2]]
            bl = length(bp)
            switch(SEX[b], {
                trans_mats = .trans_M(morh, bh, mutmat = mutmat$female)
                mm = rep(.colSums(trans_mats * bp, bl, mo_len), each = fa_len)
            }, {
                trans_pats = if (is.null(mutmat)) unlist(lapply(farh, function(fh) as.numeric(fh ==
                  bh[1, ]))) else unlist(lapply(farh, function(fh) mutmat$male[fh, bh[1, ]]))

                trans_mats = .trans_M(morh, bh[2, ], mutmat = mutmat$female)
                dim(trans_mats) = c(bl, mo_len)
                trans_mats_rep = as.numeric(do.call(rbind, rep(list(trans_mats), fa_len)))
                mm = .colSums((trans_pats * bp) * trans_mats_rep, bl, fa_len * mo_len)
            })
            likel = likel * mm
        }
        if (pivtype == 0) return(sum(likel))

        res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
            {
                pivh = dat[[c(piv, 1)]]
                pivp = dat[[c(piv, 2)]]
                pi_len = length(pivp)
                switch(SEX[piv], {
                  trans_mats = .trans_M(morh, pivh, mutmat = mutmat$female)
                  T = rep(trans_mats, each = fa_len)
                }, {
                  T = numeric(fa_len * mo_len * pi_len)
                  dim(T) = c(fa_len, mo_len, pi_len)
                  trans_mats = .trans_M(morh, pivh[2, ], mutmat = mutmat$female)
                  dim(trans_mats) = c(pi_len, mo_len)
                  for (i in seq_len(fa_len)) {
                    trans_pats = if (is.null(mutmat)) as.numeric(farh[i] == pivh[1, ]) else mutmat$male[farh[i],
                      pivh[1, ]]
                    T[i, , ] = t.default(trans_pats * trans_mats)  #TODO:make faster?
                  }
                })
                arr = as.vector(T) * as.vector(likel)
                dim(arr) = dim(T)
                res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[piv]]
                res * pivp
            })
        pivhap_update = switch(SEX[piv], dat[[c(piv, 1)]][res > 0], dat[[c(piv, 1)]][, res >
            0, drop = F])
    })
    dat[[piv]] = list(hap = pivhap_update, prob = res[res > 0])
    if (sum(res) == 0)
        attr(dat, "impossible") = TRUE
    return(dat)
}


.peel_MD <- function(dat, sub, theta, chrom, SEX) {
    far = sub[["father"]]
    mor = sub[["mother"]]
    offs = nonpiv.offs = sub[["offspringfb"]]
    piv = sub[["pivot"]]
    pivtype = sub[["pivtype"]]
    if (pivtype == 3)
        non.pivoffs = offs[offs != piv]
    farh = dat[[c(far, 1)]]
    morh = dat[[c(mor, 1)]]
    likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
    dims = dim(likel)
    fa_len = dims[1L]
    mo_len = dims[2L]

    .trans_MD <- function(parent_ph, child_haplo, theta) {
        if (theta == 0)
            return(((parent_ph[1] == child_haplo) + (parent_ph[2] == child_haplo))/2)
        parent_rec = abs(parent_ph) * sign(parent_ph[2:1])
        ((parent_ph[1] == child_haplo) * (1 - theta) + (parent_ph[2] == child_haplo) * (1 -
            theta) + (parent_rec[1] == child_haplo) * theta + (parent_rec[2] == child_haplo) *
            theta)/2
    }

    switch(chrom, AUTOSOMAL = {
        alloffs_pat = unique.default(unlist(lapply(offs, function(i) dat[[c(i, 1)]][1, ])))
        alloffs_mat = unique.default(unlist(lapply(offs, function(i) dat[[c(i, 1)]][2, ])))
        fathertrans = unlist(lapply(seq_len(fa_len), function(i) .trans_MD(farh[, i], alloffs_pat,
            theta)))
        mothertrans = unlist(lapply(seq_len(mo_len), function(j) .trans_MD(morh[, j], alloffs_mat,
            theta)))
        dim(fathertrans) = c(length(alloffs_pat), fa_len)
        dim(mothertrans) = c(length(alloffs_mat), mo_len)
        mm_init = numeric(fa_len * mo_len)
        dim(mm_init) = dims
        for (b in nonpiv.offs) {
            mm = mm_init
            bh = dat[[c(b, 1)]]
            bp = dat[[c(b, 2)]]
            b_pat = match(bh[1, ], alloffs_pat)
            b_mat = match(bh[2, ], alloffs_mat)
            for (i in seq_len(fa_len)) {
                trans_pats = fathertrans[b_pat, i]
                for (j in seq_len(mo_len)) mm[i, j] = (trans_pats * mothertrans[b_mat, j]) %*%
                  bp
            }
            likel = likel * mm
        }
        if (pivtype == 0) return(sum(likel))

        res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
            {
                pivh = dat[[c(piv, 1)]]
                pivp = dat[[c(piv, 2)]]
                pi_len = length(pivp)
                piv_pat = match(pivh[1, ], alloffs_pat)
                piv_mat = match(pivh[2, ], alloffs_mat)
                T = numeric(fa_len * mo_len * pi_len)
                dim(T) = c(fa_len, mo_len, pi_len)
                for (i in seq_len(fa_len)) {
                  trans_pats = fathertrans[piv_pat, i]
                  for (j in seq_len(mo_len)) T[i, j, ] <- trans_pats * mothertrans[piv_mat,
                    j]
                }
                arr = as.vector(T) * as.vector(likel)
                dim(arr) = dim(T)
                res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[piv]]
                res = res * pivp
            })
        pivhap_update = dat[[c(piv, 1)]][, res > 0, drop = F]
    }, X = {
        mm_init = numeric(length(likel))
        dim(mm_init) = dims
        for (b in nonpiv.offs) {
            mm = mm_init
            switch(SEX[b], for (j in seq_len(mo_len)) mm[, j] = .trans_MD(morh[, j], dat[[c(b,
                1)]], theta) %*% dat[[c(b, 2)]], {
                for (j in seq_len(mo_len)) {
                  transmother = .trans_MD(morh[, j], dat[[c(b, 1)]][2, ], theta)
                  for (i in seq_len(fa_len)) mm[i, j] = (as.numeric(farh[i] == dat[[c(b, 1)]][1,
                    ]) * transmother) %*% dat[[c(b, 2)]]
                }
            })
            likel = likel * mm
        }
        if (pivtype == 0) return(sum(likel))

        res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
            {
                pivh = dat[[c(piv, 1)]]
                pivp = dat[[c(piv, 2)]]
                pi_len = length(pivp)
                T = numeric(fa_len * mo_len * pi_len)
                dim(T) = c(fa_len, mo_len, pi_len)
                switch(SEX[piv], {
                  for (j in seq_len(mo_len)) {
                    transmother = .trans_MD(morh[, j], pivh, theta)
                    for (i in seq_len(fa_len)) T[i, j, ] = transmother
                  }
                }, {
                  for (j in seq_len(mo_len)) {
                    transmother = .trans_MD(morh[, j], pivh[2, ], theta)
                    for (i in seq_len(fa_len)) T[i, j, ] = (as.numeric(farh[i] == pivh[1, ]) *
                      transmother)
                  }
                })

                arr = as.vector(T) * as.vector(likel)
                dim(arr) = dim(T)
                res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[piv]]
                res * pivp
            })
        pivhap_update = switch(SEX[piv], dat[[c(piv, 1)]][res > 0], dat[[c(piv, 1)]][, res >
            0, drop = F])
    })
    dat[[piv]] = list(hap = pivhap_update, prob = res[res > 0])
    if (all(res == 0))
        attr(dat, "impossible") = TRUE
    return(dat)
}


.peel_MM <- function(dat, sub, theta, chrom, SEX) {
    far = sub[["father"]]
    mor = sub[["mother"]]
    offs = sub[["offspringfb"]]
    piv = sub[["pivot"]]
    pivtype = sub[["pivtype"]]
    if (pivtype == 3)
        offs = offs[offs != piv]

    farh = dat[[c(far, 1)]]
    morh = dat[[c(mor, 1)]]
    likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
    dims = dim(likel)
    fa_len = dims[1L]
    mo_len = dims[2L]
    mm_init = numeric(length(likel))
    dim(mm_init) = dims

    .trans_MM <- function(parent.haps, gamete.hap, theta) {
        # parent.haps = c(M1_1, M1_2, M2_1, M2_2)
        if (is.matrix(gamete.hap))
            vapply(seq_len(ncol(gamete.hap)), function(kol) .trans_MM(parent.haps, gamete.hap[,
                kol], theta), 1) else sum(c(all(parent.haps[c(1, 3)] == gamete.hap), all(parent.haps[c(2, 4)] == gamete.hap),
            all(parent.haps[c(1, 4)] == gamete.hap), all(parent.haps[c(2, 3)] == gamete.hap)) *
            c(1 - theta, 1 - theta, theta, theta))/2
    }

    switch(chrom, AUTOSOMAL = {
        for (b in offs) {
            mm = mm_init
            bh = dat[[c(b, 1)]]
            bp = dat[[c(b, 2)]]
            for (i in seq_len(fa_len)) {
                transfather = .trans_MM(farh[, i], bh[c(1, 3), , drop = F], theta)
                for (j in seq_len(mo_len)) mm[i, j] = (transfather * .trans_MM(morh[, j], bh[c(2,
                  4), , drop = F], theta)) %*% bp
            }
            likel <- likel * mm
        }
        if (pivtype == 0) return(sum(likel))

        res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
            {
                pivh = dat[[c(piv, 1)]]
                pivp = dat[[c(piv, 2)]]
                pi_len = length(pivp)
                T = numeric(fa_len * mo_len * pi_len)
                dim(T) = c(fa_len, mo_len, pi_len)
                for (i in seq_len(fa_len)) {
                  transfather = .trans_MM(farh[, i], pivh[c(1, 3), , drop = F], theta)
                  for (j in seq_len(mo_len)) T[i, j, ] <- transfather * .trans_MM(morh[, j],
                    pivh[c(2, 4), , drop = F], theta)
                }
                arr = as.vector(T) * as.vector(likel)
                dim(arr) = dim(T)
                res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[piv]]
                res = res * pivp
            })
        pivhap_update = dat[[c(piv, 1)]][, res > 0, drop = F]
    }, X = {
        for (b in offs) {
            mm = mm_init
            bh = dat[[c(b, 1)]]
            bp = dat[[c(b, 2)]]
            if (SEX[b] == 1) for (j in seq_len(mo_len)) mm[, j] = .trans_MM(morh[, j], bh,
                theta) %*% bp else for (j in seq_len(mo_len)) {
                transmother = .trans_MM(morh[, j], bh[c(2, 4), , drop = F], theta)
                for (i in seq_len(fa_len)) mm[i, j] = (as.numeric(farh[1, i] == bh[1, ] & farh[2,
                  i] == bh[3, ]) * transmother) %*% bp
            }
            likel = likel * mm
        }
        if (pivtype == 0) return(sum(likel))

        res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
            {
                pivh = dat[[c(piv, 1)]]
                pivp = dat[[c(piv, 2)]]
                pi_len = length(pivp)
                T = numeric(fa_len * mo_len * pi_len)
                dim(T) = c(fa_len, mo_len, pi_len)
                if (SEX[piv] == 1) {
                  for (j in seq_len(mo_len)) {
                    transmother = .trans_MM(morh[, j], pivh, theta)
                    for (i in seq_len(fa_len)) T[i, j, ] = transmother
                  }
                } else {
                  for (j in seq_len(mo_len)) {
                    transmother = .trans_MM(morh[, j], pivh[c(2, 4), ], theta)
                    for (i in seq_len(fa_len)) T[i, j, ] = (as.numeric(farh[1, i] == pivh[1,
                      ] & farh[2, i] == pivh[3, ]) * transmother)
                  }
                }
                arr = as.vector(T) * as.vector(likel)
                dim(arr) = dim(T)
                res = .colSums(arr, fa_len * mo_len, pi_len)  #sum for each entry of haps[[piv]]
                res * pivp
            })
        pivhap_update = dat[[c(piv, 1)]][, res > 0, drop = F]
    })
    dat[[piv]] = list(hap = pivhap_update, prob = res[res > 0])
    if (sum(res) == 0)
        attr(dat, "impossible") = TRUE
    return(dat)
}



###### OTHER AUXILIARY FUNCTIONS

.informative = function(x, marker, marker2 = NULL) {
    # Trim pedigree by removing leaves without genotypes, and also remove completely
    # uninformative subnucs.
    if (!is.null(marker2))
        marker = marker + marker2
    if (all(marker[, 1] > 0))
        return(list(subnucs = x$subnucs, newfounders = numeric(0)))
    newfounders = numeric(0)
    new_subnucs = list()
    p = x$pedigree
    is_miss = marker[, 1] == 0 & marker[, 2] == 0
    is_miss[loop_int <- .internalID(x, x$loop_breakers)] = F  # works (and quick) also if no loops.
    is_uninf_leaf = is_miss & !p[, "ID"] %in% p[, c("FID", "MID")]
    is_uninf_fou = is_miss & seq_len(x$nInd) %in% x$founders

    for (sub in x$subnucs) {
        fa = sub[["father"]]
        mo = sub[["mother"]]
        offs = sub[["offspringfb"]]
        pivot = sub[["pivot"]]
        sub[["offspringfb"]] = offs[!is_uninf_leaf[offs]]
        noffs = length(sub[["offspringfb"]])
        switch(sub[["pivtype"]], {
            if (noffs == 0 && is_uninf_fou[mo]) sub = NULL
        }, {
            if (noffs == 0 && is_uninf_fou[fa]) sub = NULL
        }, {
            if (noffs == 1 && is_uninf_fou[fa] && is_uninf_fou[mo]) {
                newfounders = c(newfounders, pivot)
                sub = NULL
            }
        })
        if (!is.null(sub)) {
            new_subnucs = c(new_subnucs, list(sub))
            is_uninf_fou[pivot] = FALSE  #added in v0.8-1 to correct a bug marking certain 'middle' subnucs uninformative
        }
    }
    list(subnucs = new_subnucs, newfounders = newfounders)
}

