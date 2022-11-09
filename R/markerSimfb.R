#' Marker simulation
#'
#' Simulates marker genotypes conditional on the pedigree structure, affection
#' statuses and disease model.
#'
#' This implements (with various time savers) the algorithm used in SLINK of the
#' LINKAGE/FASTLINK suite. If \code{partialmarker} is NULL, genotypes are
#' simulated by simple gene dropping, using \code{simpleSim}.
#'
#' @param x a \code{linkdat} object
#' @param N a positive integer: the number of markers to be simulated
#' @param available a vector containing IDs of the available individuals, i.e.
#'   those whose genotypes should be simulated. By default, all individuals are
#'   included.
#' @param alleles a vector containing the alleles for the marker to be
#'   simulation. If a single integer is given, this is interpreted as the number
#'   of alleles, and the actual alleles as \code{1:alleles}. Must be NULL if
#'   \code{partialmarker} is non-NULL.
#' @param afreq a vector of length 2 containing the population frequencies for
#'   the marker alleles. Must be NULL if \code{partialmarker} is non-NULL.
#' @param partialmarker Either NULL (resulting in unconditional simulation), a
#'   marker object (on which the simulation should be conditioned) or the index
#'   of an existing marker of \code{x}.
#' @param loop_breakers a numeric containing IDs of individuals to be used as
#'   loop breakers. Relevant only if the pedigree has loops, and only if
#'   \code{partialmarker} is non-NULL. See \code{breakLoopsfb}.
#' @param eliminate A non-negative integer, indicating the number of iterations
#'   in the internal genotype-compatibility algorithm. Positive values can save
#'   time if \code{partialmarker} is non-NULL and the number of alleles is
#'   large.
#' @param seed NULL, or a numeric seed for the random number generator.
#' @param verbose a logical.
#' @return a \code{linkdat} object equal to \code{x} except its
#'   \code{markerdata} entry, which consists of the \code{N} simulated markers.
#' @seealso \code{simpleSim}, \code{linkage.power}
#' @references G. M. Lathrop, J.-M. Lalouel, C. Julier, and J. Ott,
#'   \emph{Strategies for Multilocus Analysis in Humans}, PNAS 81(1984), pp.
#'   3443-3446.
#' @import utils
#'
#' @export
markerSimfb <- function(x, N = 1, available = NULL, alleles = NULL, afreq = NULL, partialmarker = NULL,
    loop_breakers = NULL, eliminate = 0, seed = NULL, verbose = TRUE) {

    if (!is.linkdat(x) && !is.linkdat.list(x))
        stop("x must be either a 'linkdat' object, a 'singletonfb' object, or a list of such")

    if (!is.null(seed))
        set.seed(seed)

    # if input is a list of linkdat objects: Apply markerSimfb recursively
    if (is.linkdat.list(x))
        return(lapply(x, function(xi) markerSimfb(xi, N = N, available = intersect(xi$orig.ids, available),
            alleles = alleles, afreq = afreq, partialmarker = partialmarker,
            loop_breakers = loop_breakers, eliminate = eliminate, verbose = verbose)))

    starttime = proc.time()

    # Reorder if necessary
    if(!.check_parentsfb_before_children(x)) {
        if(verbose) cat("Note: Changing the internal order so that all parentsfb precede their children.\n\n")
        x = .reorder_parentsfb_before_children(x)
    }

    likel_counter = 0
    assert_that(.is.natural(N))
    if (!is.null(x$loop_breakers))
        stop("Linkdat objects with pre-broken loops are not allowed as input to the `markerSimfb` function.")
    if(is.null(available)) available = x$orig.ids

    if (!is.null(partialmarker)) {
        if (!is.null(alleles) || !is.null(afreq))
            stop("When 'partialmarker' is non-NULL, both 'alleles' and 'afreq' must be NULL.")

        if (inherits(partialmarker, "marker")) {
            if (nrow(partialmarker) != x$nInd)
                stop("Partial marker does not fit the pedigree.")
        } else if (.is.natural(partialmarker) && partialmarker <= x$nMark) {
            partialmarker = x$markerdata[[partialmarker]]
        } else stop("The 'partialmarker' must be a 'marker' object, or a single integer indicating an existing marker of 'x'.")

        if (is.null(attr(partialmarker, "mutmat"))) {
            err = mendelianCheckfb(SetMarkersfb(x, partialmarker), verbose = F)
            if (length(err) > 0)
                stop("Mendelian error in the given partial marker.")
        }
    } else {
        if (is.null(alleles))
            stop("Please specify marker alleles.")
        if (is.numeric(alleles) && length(alleles) == 1)
            alleles = seq_len(alleles)
        partialmarker = markerfb(x, alleles = alleles, afreq = afreq)
    }
    m = partialmarker
    alleles = attr(m, "alleles")
    afreq = attr(m, "afreq")
    mutmat = attr(m, "mutmat")
    chromattr = attr(m, "chrom")
    chrom = if (identical(23L, as.integer(chromattr)))
        "X" else "AUTOSOMAL"

    if (all(m == 0))
        return(simpleSimfb(x, N, alleles = alleles, afreq = afreq, available = available, Xchrom = (chrom ==
            "X"), mutmat = mutmat, seed = seed, verbose = verbose))

    allgenos = allGenotypes(nall <- attr(m, "nalleles"))
    mutations = !is.null(mutmat)

    gridlist = geno.grid.subset(x, m, x$orig.ids, chrom, make.grid = F)

    # Forced genotypes:
    forcedTF = (m[, 1] == 0 | m[, 2] == 0) & (lengths(gridlist) == 1)
    m.unforced = m  # saving a copy for use in verbose output below
    for (id in (1:x$nInd)[forcedTF]) m[id, ] = allgenos[gridlist[[id]], ]

    if (verbose) {
        cat(ifelse(chrom == "AUTOSOMAL", "Autosomal", "X-linked"), "marker locus.\n")
        cat(sprintf("Simulating genotypes for individual%s %s.\n", ifelse(length(available) ==
            1, "", "s"), .prettycat(available, "and")))
        cat("\nAlleles and frequencies:\n")
        print(structure(round(afreq, 3), names = alleles))
        if (mutations) {
            cat("\nMutation matrices:\n")
            print(mutmat)
        } else cat("\n")
        cat("Conditioning on the following genotypes:\n")
        print(data.frame(ID = x$orig.ids, GENO = .prettyMarkers(m.unforced, missing = "-",
            singleCol = TRUE, sep = "/", sex = x$pedigree[, "SEX"])))
        if (any(forcedTF)) {
            cat("\nForced genotypes:\n")
            for (id in (1:x$nInd)[forcedTF]) {
                allelchars = alleles[m[id, ]]
                if (chrom == "X")
                  allelchars = allelchars[1]
                cat(sprintf("Individual %d: %s\n", x$orig.ids[id], paste(allelchars, collapse = "/")))
            }
        }
    }

    # Making copies of x and m before possible loop breaking (used to determine simulation
    # strategy)
    xorig = SetMarkersfb(x, NULL)
    morig = m

    if (loops <- x$hasLoops) {
        orig_ids = x$orig.ids
        x = breakLoopsfb(SetMarkersfb(x, m), loop_breakers = loop_breakers, verbose = verbose)
        m = x$markerdata[[1]]
        loop_breakers = x$loop_breakers[, 1]
        gridlist = gridlist[sort.int(match(c(orig_ids, loop_breakers), orig_ids))]
    }

    ngrid = lengths(gridlist)
    ped = x$pedigree
    SEX = ped[, "SEX"]

    ##### Determine simulation strategy #### Note: Using original x and m in this section (i.e.
    ##### before loop breaking)

    # Individuals that are typed (or forced - see above). Simulations condition on these.
    typedTF = (morig[, 1] != 0 | morig[, 2] != 0)
    typed = xorig$orig.ids[typedTF]

    # Target individuals: untyped individuals that we seek to simulate
    targets = .mysetdiff(available, typed)
    untyped_breakers = if (loops)
        .mysetdiff(loop_breakers, typed) else NULL

    # Method 2: Compute joint dist of some target individuals, brute force on the remaining
    hardsim.method2 = unique.default(c(untyped_breakers, targets))
    hardsim.method2_int = .internalID(x, hardsim.method2)
    method2 = .optimal.precomputation(hardsim.method2_int, N, gridlist, chrom, SEX = SEX)

    ### Method 3: Extend target to ancestorsfb of targets.
    targets.plus.ancestorsfb = .mysetdiff(c(targets, ancestorsfb(xorig, id = targets)), typed)

    # Only ancestorsfb of typed individuals are hard; the others can be simple-dropped
    ancestorsfb.of.typed = ancestorsfb(xorig, id = typed)

    hardsim.method3 = intersect(targets.plus.ancestorsfb, ancestorsfb.of.typed)
    hardsim.method3 = unique.default(c(untyped_breakers, hardsim.method3))

    hardsim.method3_int = .internalID(x, hardsim.method3)
    method3 = .optimal.precomputation(hardsim.method3_int, N, gridlist, chrom, SEX = SEX)

    if (method2$calls <= method3$calls) {
        joint_int = method2$id_int
        bruteforce_int = .mysetdiff(hardsim.method2_int, joint_int)
        simpledrop = numeric()
    } else {
        joint_int = method3$id_int
        bruteforce_int = .mysetdiff(hardsim.method3_int, joint_int)
        simpledrop = .mysetdiff(targets.plus.ancestorsfb, ancestorsfb.of.typed)
    }

    simpledrop_int = .internalID(x, simpledrop)
    simple.founders_int = intersect(simpledrop_int, x$founders)
    simple.nonfounders_int = intersect(simpledrop_int, x$nonfounders)

    if (length(simple.nonfounders_int) > 0) {
        # Ensure sensible ordering of nonfounders (for gene dropping)
        done = c(.internalID(x, typed), joint_int, bruteforce_int, simple.founders_int)
        v = simple.nonfounders_int
        v.ordered = numeric()
        while (length(v) > 0) {
            i = match(T, (ped[v, "FID"] %in% done) & (ped[v, "MID"] %in% done))
            if (is.na(i))
                stop("Could not determine sensible order for gene dropping. Bug report to magnusdv@medisin.uio.no is appreciated!")
            v.ordered = c(v.ordered, v[i])
            done = c(done, v[i])
            v = v[-i]
        }
        simple.nonfounders_int = v.ordered
    }

    if (verbose) {
        cat("\nSimulation strategy:\n")

        .tostring = function(v) if (length(v) > 0)
            .prettycat(x$orig.ids[v], "and") else "None"
        cat(sprintf("Pre-computed joint distribution: %s.\n", .tostring(joint_int)))
        cat(sprintf("Brute force conditional simulation: %s.\n", .tostring(bruteforce_int)))
        cat(sprintf("Hardy-Weinberg sampling (founders): %s.\n", .tostring(simple.founders_int)))
        cat(sprintf("Simple gene dropping: %s.\n", .tostring(simple.nonfounders_int)))
        cat(sprintf("Required likelihood computations: %d\n", min(method2$calls, method3$calls)))
    }
    # create initial marker matrix: two columns per marker
    markers = rep.int(m, N)
    dim(markers) = c(x$nInd, 2 * N)
    odd = seq_len(N) * 2 - 1

    if (length(joint_int) > 0) {
        allgenos_row_grid = t.default(fast.grid(gridlist[joint_int]))  #Cartesian product. Each row contains 'init' row numbers of allgenos.
        jointp = apply(allgenos_row_grid, 2, function(rownrs) {
            partial = m
            partial[joint_int, ] = allgenos[rownrs, ]
            likelihood.linkdat(x, locus1 = partial, eliminate = eliminate)
        })
        likel_counter = likel_counter + length(jointp)
        if (identical(sum(jointp), 0))
            stop("When trying to pre-compute joint probabilities: All probabilities zero. Mendelian error?")

        # fill the rows of the 'joint' individuals
        sample_rows = allgenos_row_grid[, suppressWarnings(sample.int(length(jointp), size = N,
            replace = TRUE, prob = jointp))]
        markers[joint_int, odd] = allgenos[sample_rows, 1]
        markers[joint_int, odd + 1] = allgenos[sample_rows, 2]
    }

    if (length(bruteforce_int) > 0) {
        for (i in bruteforce_int) {
            gridi = gridlist[[i]]
            rowsample = unlist(lapply(2 * seq_len(N), function(mi) {
                partial = m
                partial[] = markers[, c(mi - 1, mi)]  # preserves all attributes of the m.
                probs = unlist(lapply(gridi, function(r) {
                  partial[i, ] = allgenos[r, ]
                  likelihood.linkdat(x, locus1 = partial, eliminate = eliminate)
                }))
                if (sum(probs) == 0) {
                  print(cbind(ped, partial))
                  stop("\nIndividual ", x$orig.ids[i], ": All genotype probabilities zero. Mendelian error?")
                }
                sample(gridi, size = 1, prob = probs)
            }))
            markers[i, odd] = allgenos[rowsample, 1]
            markers[i, odd + 1] = allgenos[rowsample, 2]
        }
        likel_counter = likel_counter + N * sum(ngrid[bruteforce_int])
    }

    if (length(simpledrop) > 0) {
        loopbr_int = .internalID(x, x$loop_breakers[, 1])  #integer(0) if no loops
        loopbr_dup_int = .internalID(x, x$loop_breakers[, 2])

        if (chrom == "AUTOSOMAL")
            markers[simple.founders_int, ] = sample.int(nall, size = 2 * N * length(simple.founders_int),
                replace = TRUE, prob = afreq) else for (f in simple.founders_int) markers[f, ] = switch(SEX[f], rep(sample.int(nall,
            size = N, replace = TRUE, prob = afreq), each = 2), sample.int(nall, size = 2 *
            N, replace = TRUE, prob = afreq))

        markers[loopbr_dup_int, ] = markers[loopbr_int, ]  # Genotypes of the duplicated individuals. Some of these may be ungenotyped...save time by excluding these?

        for (id in simple.nonfounders_int) {
            fa = ped[id, "FID"]
            mo = ped[id, "MID"]
            if (chrom == "AUTOSOMAL") {
                paternal = markers[fa, odd + .rand01(N)]
                maternal = markers[mo, odd + .rand01(N)]
                if (mutations) {
                  paternal = vapply(paternal, function(a) sample.int(nall, size = 1, prob = mutmat$male[a,
                    ]), 1)
                  maternal = vapply(maternal, function(a) sample.int(nall, size = 1, prob = mutmat$female[a,
                    ]), 1)
                }
            } else {
                maternal = markers[mo, odd + .rand01(N)]
                if (mutations)
                  maternal = vapply(maternal, function(a) sample.int(nall, size = 1, prob = mutmat$female[a,
                    ]), 1)

                if (SEX[id] == 1)
                  paternal = maternal  # if boy, only maternal
                else {
                  paternal = markers[fa, odd]  # if girl, fathers allele is forced
                  if (mutations)
                    paternal = vapply(paternal, function(a) sample.int(nall, size = 1, prob = mutmat$male[a,
                      ]), 1)
                }
            }
            markers[id, odd] = paternal
            markers[id, odd + 1] = maternal
        }
    }

    if (loops) {
        markers = markers[-match(x$loop_breakers[, 2], x$orig.ids), ]
        x = tieLoopsfb(x)
    }

    # removing genotypes for individuals that are i) originally untyped and ii) unavailable
    typedTF[forcedTF] = F
    unavailable = !(x$orig.ids %in% available)
    markers[!typedTF & unavailable, ] = 0
    attrib = attributes(partialmarker)
    attrib$name = NA
    markerdata_list = lapply(seq_len(N), function(k) {
        mk = markers[, c(2 * k - 1, 2 * k)]
        attributes(mk) = attrib
        mk
    })
    class(markerdata_list) = "markerdata"
    x = SetMarkersfb(x, markerdata_list)
    seconds = (proc.time() - starttime)[["elapsed"]]
    if (verbose)
        cat(sprintf("\n%d markers simulated.\nNumber of calls to the likelihood function: %d.\nTotal time used: %f seconds.\n",
            x$nMark, likel_counter, seconds))
    x
}

.optimal.precomputation = function(target_int, Nsim, gridlist, chrom, SEX = NULL) {
    if (length(target_int) == 0)
        return(list(calls = 0, id_int = target_int))
    if (chrom == "AUTOSOMAL") {
        T = length(target_int)
        ngrid_target = lengths(gridlist[target_int])
        callsCum = sapply(1:T, function(ci) prod(ngrid_target[seq_len(ci)]) + Nsim * sum(ngrid_target[seq_len(T -
            ci) + ci]))
        minimum_index = which.min(callsCum)
        opt = list(calls = callsCum[minimum_index], id_int = target_int[seq_len(minimum_index)])
    } else if (chrom == "X") {
        males = target_int[SEX[target_int] == 1]
        females = target_int[SEX[target_int] == 2]
        ngrid_m = lengths(gridlist[males], use.names = F)
        ngrid_f = lengths(gridlist[females], use.names = F)
        M = length(males)
        F = length(females)
        # Find optimal 'init' values for males/females (fewest likelihood calls)
        callsCum = matrix(nrow = M + 1, ncol = F + 1)
        for (ma in 0:M) for (fe in 0:F) callsCum[ma + 1, fe + 1] = prod(ngrid_m[seq_len(ma)]) *
            prod(ngrid_f[seq_len(fe)]) + Nsim * sum(c(ngrid_m[seq_len(M - ma) + ma], ngrid_f[seq_len(F -
            fe) + fe]))

        minimum_index = arrayInd(which.min(callsCum), dim(callsCum))
        id_int = c(males[seq_len(minimum_index[1] - 1)], females[seq_len(minimum_index[2] -
            1)])
        opt = list(calls = callsCum[minimum_index], id_int = id_int)
    }
    return(opt)
}






#' Unconditional marker simulation
#'
#' Unconditional simulation of unlinked markers
#'
#' This simulation is done by distributing alleles randomly to all founders,
#' followed by unconditional gene dropping down throughout the pedigree (i.e.
#' for each non-founder a random allele is selected from each of the parentsfb).
#' Finally the genotypes of any individuals not included in \code{available} are
#' removed.
#'
#' @param x a \code{linkdat} object
#' @param N a positive integer: the number of markers to be simulated
#' @param alleles a vector containing the allele names. If missing, the alleles
#'   are taken to be \code{seq_along(afreq)}.
#' @param afreq a vector of length 2 containing the population frequencies for
#'   the alleles. If missing, the alleles are assumed equifrequent.
#' @param available a vector containing IDs of the available individuals, i.e.
#'   those whose genotypes should be simulated.
#' @param Xchrom a logical: X linked markers or not?
#' @param mutmat a mutation matrix, or a list of two such matrices named
#'   'female' and 'male'. The matrix/matrices must be square, with the allele
#'   labels as dimnames, and each row must sum to 1 (after rounding to 3
#'   decimals).
#' @param seed NULL, or a numeric seed for the random number generator.
#' @param verbose a logical.
#' @return a \code{dat} object equal to \code{x} in all respects except its
#'   \code{markerdata} entry, which consists of the \code{N} simulated markers.
#' @seealso \code{markerSimfb}, \code{linkageSim}
#'
#'
#' @export
simpleSimfb = function(x, N, alleles, afreq, available, Xchrom = FALSE, mutmat = NULL, seed = NULL,
    verbose = T) {
    starttime = proc.time()
    if (missing(alleles)) {
        if (missing(afreq))
            stop("Both 'alleles' and 'afreq' cannot be missing")
        alleles = seq_along(afreq)
    }
    nall = length(alleles)
    if (missing(afreq))
        afreq = rep(1, nall)/nall
    if (variableSNPfreqs <- (nall == 2 && length(afreq) != 2 && !Xchrom))
        afreq = rep(afreq, length = N)
    if (missing(available))
        available = x$orig.ids
    mutations = !is.null(mutmat)
    if (mutations) {
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

    # Reorder if necessary
    if(!.check_parentsfb_before_children(x)) {
        if(verbose) cat("Note: Changing the internal order so that all parentsfb precede their children.\n\n")
        x = .reorder_parentsfb_before_children(x)
    }

    if (verbose) {
        plural.s = ifelse(length(available) == 1, "", "s")
        cat(ifelse(!Xchrom, "Autosomal", "X-linked"), "marker locus.\n")
        cat(sprintf("Unconditional simulating of genotypes for individual%s %s.\n", plural.s,
            .prettycat(available, "and")))
        cat("\nAlleles and frequencies:\n")
        if (variableSNPfreqs)
            cat("SNPs with allele 1 frequencies", paste(utils::head(afreq, 5), collapse = ", "), ifelse(N >
                5, "...\n", "\n")) else print(structure(round(afreq, 3), names = alleles))
        if (mutations) {
            cat("\nMutation matrices:\n")
            print(mutmat)
        } else cat("\n")
    }

    ped = x$pedigree
    m = matrix(0, ncol = 2 * N, nrow = x$nInd)
    odd = seq_len(N) * 2 - 1

    if (!is.null(seed))
        set.seed(seed)
    if (Xchrom) {
        for (f in x$founders) {
            if (ped[f, "SEX"] == 1)
                m[f, ] = rep(sample.int(nall, size = N, replace = TRUE, prob = afreq), each = 2)
            if (ped[f, "SEX"] == 2)
                m[f, ] = sample.int(nall, size = 2 * N, replace = TRUE, prob = afreq)
        }
        for (id in x$nonfounders) {
            fa = ped[id, "FID"]
            mo = ped[id, "MID"]
            maternal = m[mo, odd + .rand01(N)]
            if (mutations) {
                maternal = vapply(maternal, function(a) sample.int(nall, size = 1, prob = mutmat$female[a, ]), 1)
            }
            if (ped[id, "SEX"] == 1) {
                paternal = maternal  # if boy, only maternal
            } else {
                paternal = m[fa, odd]  # if girl, fathers allele is forced
                if (mutations)
                  paternal = vapply(paternal, function(a) sample.int(nall, size = 1, prob = mutmat$male[a,
                    ]), 1)
            }
            m[id, odd] = paternal
            m[id, odd + 1] = maternal
        }
    } else {
        size = 2 * length(x$founders)
        allelsamp = if (variableSNPfreqs)
            unlist(lapply(afreq, function(f) sample.int(2, size, replace = TRUE, prob = c(f,
                1 - f)))) else sample.int(nall, size = N * size, replace = TRUE, prob = afreq)
        m[x$founders, ] = allelsamp
        for (id in x$nonfounders) {
            fa = ped[id, "FID"]
            mo = ped[id, "MID"]
            paternal = m[fa, odd + .rand01(N)]
            maternal = m[mo, odd + .rand01(N)]
            if (!is.null(mutmat)) {
                paternal = vapply(paternal, function(a) sample.int(nall, size = 1, prob = mutmat$male[a,
                  ]), 1)
                maternal = vapply(maternal, function(a) sample.int(nall, size = 1, prob = mutmat$female[a,
                  ]), 1)
            }
            m[id, odd] = paternal
            m[id, odd + 1] = maternal
        }
    }

    m[!x$orig.ids %in% available, ] = 0
    if (variableSNPfreqs) {
        attrib = attributes(markerfb(x, alleles = alleles, afreq = NULL, chrom = NA, mutmat = mutmat,
            missing = 0))
        frqs = as.vector(rbind(afreq, 1 - afreq))
        markerdata_list = lapply(odd, function(k) {
            mk = m[, c(k, k + 1)]
            atr = attrib
            atr$afreq = frqs[c(k, k + 1)]
            attributes(mk) = atr
            mk
        })
    } else {
        attrib = attributes(markerfb(x, alleles = alleles, afreq = afreq, chrom = ifelse(Xchrom,
            23, NA), mutmat = mutmat, missing = 0))
        markerdata_list = lapply(odd, function(k) {
            mk = m[, c(k, k + 1)]
            attributes(mk) = attrib
            mk
        })
    }
    x = SetMarkersfb(x, structure(markerdata_list, class = "markerdata"))
    if (verbose) {
        seconds = (proc.time() - starttime)[["elapsed"]]
        cat(sprintf("%d markers simulated.\nNumber of calls to the likelihood function: 0.\nTotal time used: %f seconds.\n",
            x$nMark, seconds))
    }
    x
}

