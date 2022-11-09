#' Set, change or display the model parameters for 'linkdat' objects
#'
#' Functions to set, change and display model parameters involved in parametric
#' linkage analysis.
#'
#' @param x in \code{setModel}: a \code{\link{linkdat}} object. In
#'   \code{print.linkdat.model}: a \code{linkdat.model} object.
#' @param model NULL, or an object of class \code{linkdat.model}, namely a list
#'   with elements \code{chrom}, \code{penetrances} and \code{dfreq}.  In the
#'   \code{setModel} function, the \code{model} argument can be one of the
#'   integers 1-4, with the following meanings:
#'
#'   1 = autosomal dominant; fully penetrant, dfreq=1e-5
#'
#'   2 = autosomal recessive; fully penetrant, dfreq=1e-5
#'
#'   3 = X-linked dominant; fully penetrant, dfreq=1e-5
#'
#'   4 = X-linked recessive; fully penetrant, dfreq=1e-5
#' @param chrom a character, either 'AUTOSOMAL' or 'X'. Lower case versions are
#'   allowed and will be converted automatically.
#' @param penetrances if \code{chrom=='AUTOSOMAL'}: a numeric of length 3 -
#'   \code{(f0, f1, f2)} - where \code{fi} is the probability of being affected
#'   given \code{i} disease alleles.
#'
#'   If \code{chrom=='X'}: a list of two vectors, containing the penetrances for
#'   each sex: \code{penetrances = list(male=c(f0, f1), female=c(f0, f1, f2))}.
#' @param dfreq the population frequency of the disease allele.
#' @param ... further parameters
#' @return \code{setModel} returns a new \code{linkdat} object, whose
#'   \code{model} entry is a \code{linkdat.model} object: A list containing the
#'   given \code{chrom}, \code{penetrances} and \code{dfreq}.
#' @seealso \code{\link{linkdat}}
#'
#'
#' @export
setModel = function(x, model = NULL, chrom = NULL, penetrances = NULL, dfreq = NULL) {
    assert_that(is.linkdat(x))
    if (!is.null(chrom))
        stopifnot(is.character(chrom))
    if (!is.null(penetrances))
        stopifnot(is.numeric(unlist(penetrances)), max(unlist(penetrances)) <= 1, min(unlist(penetrances)) >=
            0)
    if (!is.null(dfreq))
        stopifnot(is.numeric(dfreq), length(dfreq) == 1, dfreq >= 0, dfreq <= 1)

    if (is.numeric(model)) {
        stopifnot(model %in% 1:4)
        model = switch(model, list(chrom = "AUTOSOMAL", penetrances = c(0, 1, 1), dfreq = 1e-05),
            list(chrom = "AUTOSOMAL", penetrances = c(0, 0, 1), dfreq = 1e-05), list(chrom = "X",
                penetrances = list(male = c(0, 1), female = c(0, 1, 1)), dfreq = 1e-05), list(chrom = "X",
                penetrances = list(male = c(0, 1), female = c(0, 0, 1)), dfreq = 1e-05))
    }
    if (is.null(model) && !is.null(x$model))
        # If no model is given, but x already has one, use this as template
    model = x$model
    hasmodel = !is.null(model)

    if (is.null(chrom))
        chrom = ifelse(hasmodel, model$chrom, "AUTOSOMAL") else chrom <- match.arg(toupper(chrom), c("AUTOSOMAL", "X"))
    if (is.null(dfreq))
        dfreq = ifelse(hasmodel, model$dfreq, 1e-05)
    if (is.null(penetrances))
        if (hasmodel)
            penetrances = model$penetrances else stop("No penetrance values given.") else {
        switch(chrom, AUTOSOMAL = {
            if (is.character(penetrances)) {
                mod <- match.arg(tolower(penetrances), c("dominant", "recessive"))
                penetrances <- switch(mod, dominant = c(0, 1, 1), recessive = c(0, 0, 1))
            }
            if (length(penetrances) != 3) stop("For autosomal models, the penetrance parameter must be a vector of the form: c(f_0, f_1, f_2).")
            names(penetrances) <- c("f0", "f1", "f2")
        }, X = {
            if (is.character(penetrances)) {
                mod <- match.arg(tolower(penetrances), c("dominant", "recessive"))
                penetrances <- switch(mod, dominant = list(c(0, 1), c(0, 1, 1)), recessive = list(c(0,
                  1), c(0, 0, 1)))
            }
            if (any(length(penetrances) != 2, length(penetrances[[1]]) != 2, length(penetrances[[2]]) !=
                3)) stop("For X-linked models, the penetrance parameter must be a list of the form: list(c(f0_m, f1_m), c(f0_f, f1_f, f2_f)).")
            names(penetrances) <- c("male", "female")
            names(penetrances$male) <- c("f0_m", "f1_m")
            names(penetrances$female) <- c("f0_f", "f1_f", "f2_f")
        })
    }

    # Collect all model information
    x$model = structure(list(chrom = chrom, penetrances = penetrances, dfreq = dfreq), class = "linkdat.model")
    return(invisible(x))
}

#' @export
#' @rdname setModel
print.linkdat.model <- function(x, ...) {
    model = x
    switch(model$chrom, AUTOSOMAL = cat("Autosomal inheritance with penetrances: (f0, f1, f2) =",
        paste("(", paste(model$penetrances, collapse = ", "), ")", sep = ""), "\n"), X = cat("X-linked inheritance with penetrances:\n\tMales: (f0, f1) =",
        paste("(", paste(model$penetrances$male, collapse = ", "), ")", sep = ""), "\n\tFemales: (f0, f1, f2) =",
        paste("(", paste(model$penetrances$female, collapse = ", "), ")", sep = ""), "\n"))
    cat("Disease allele frequency:", model$dfreq, "\n")
}
