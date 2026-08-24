<img src="genis.png" align="left" width="80">

# fbnet

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/fbnet)](https://CRAN.R-project.org/package=fbnet)
[![](https://cranlogs.r-pkg.org/badges/grand-total/fbnet?color=blue)](https://cran.r-project.org/package=fbnet)

<!-- badges: end -->

fbnet computes likelihood ratios for kinship testing and human identification by
representing a pedigree as a Bayesian network. Every allele of every individual
becomes a node, edges encode Mendelian transmission, and queries are answered by
variable elimination rather than by peeling the pedigree directly. That
representation is what makes the package suited to large-size queries against
databases of missing individuals: the elimination order is computed once and
reused across queries, and unobserved parts of the network can be pruned before
the expensive step.

The package implements the core kinship engine of *GENis*, the open-source
forensic DNA database system developed by Fundación Sadosky
(https://github.com/fundacion-sadosky/genis).

Maintenance status: this repository is kept for the CRAN release. Active
development of new releases happens at https://github.com/chernolabs/fbnet-cran.

## Installation

``` r
install.packages("fbnet")
```

The development version:

``` r
# install.packages("remotes")
remotes::install_github("MarsicoFL/fbnet")
```

Pedigrees are built with `pedtools`, and `forrel` is useful for simulating
profiles:

``` r
install.packages(c("pedtools", "forrel"))
```

## A first example

A mother and a child, two markers with four alleles each. The mother is
genotyped `3/4` at both markers; the child is the query person, and we want the
probability distribution of the alleles they carry.

``` r
library(fbnet)
library(pedtools)

ped <- nuclearPed(1)
m1  <- marker(ped, afreq = c("1" = 0.1, "2" = 0.2, "3" = 0.3, "4" = 0.4), name = "M1")
m2  <- marker(ped, afreq = c("1" = 0.1, "2" = 0.2, "3" = 0.3, "4" = 0.4), name = "M2")
ped <- setMarkers(ped, list(m1, m2))
ped <- setGenotype(ped, id = "2", marker = 1, geno = "3/4")
ped <- setGenotype(ped, id = "2", marker = 2, geno = "3/4")

pbn  <- from_pedtools(ped, QP = "3")
bnet <- buildBN(pbn, QP = "3")
bn1  <- buildCPTs(bnet)
resQ <- velim.bn(bn1, ordMethod = "min_fill", verbose = FALSE)
```

`resQ` holds one table per queried allele, named `<individual>_<marker>_<p|m>`,
where `p` is the paternal and `m` the maternal allele. The maternal allele of the
child is constrained by the evidence, since the mother can only have transmitted
a 3 or a 4:

``` r
resQ[["3_M1_m"]]
#>   3_M1_m prob
#> 1      3 0.12
#> 2      4 0.12
```

The paternal allele has no evidence bearing on it, so its distribution is the
population one:

``` r
resQ[["3_M1_p"]]
#>   3_M1_p prob
#> 1      1  0.1
#> 2      2  0.2
#> 3      3  0.3
#> 4      4  0.4
```

Reading those two tables side by side is the point of the method: the network
reports where the evidence acts and where it does not, one allele at a time,
instead of collapsing everything into a single number.

## Likelihood ratios

Given an observed genotype for the query person, `reportLR` returns the LR per
marker and for the whole profile. The genotypes go in as a character matrix,
with individuals in rows and markers in columns.

``` r
geno <- matrix(c("3/4", "3/4"), nrow = 1,
               dimnames = list("3", c("M1", "M2")))

reportLR(bnet, resQ, geno = geno)
#> $systemLR
#>        L0    L1       LR
#> M1 0.0576 0.084 1.458333
#> M2 0.0576 0.084 1.458333
#>
#> $LR
#> [1] 2.126736
```

The per-marker value can be checked by hand. Under the pedigree hypothesis the
child inherits a 3 or a 4 from the mother with probability 1/2 each and the
second allele from the population, so P(3/4) = 0.5 x 0.4 + 0.5 x 0.3 = 0.35.
Under the alternative the child is an unrelated individual, so
P(3/4) = 2 x 0.3 x 0.4 = 0.24. The ratio is 0.35 / 0.24 = 1.4583, and the profile
LR is that value squared, the two markers being independent.

The denominator on its own, the random match probability of the queried
genotype, is available directly:

``` r
getQSetRMP(bnet, list("3_M1" = c("3", "4"), "3_M2" = c("3", "4")))
#> [1] 0.0576
```

## Pedigrees from pedtools and forrel

`from_pedtools()` accepts any `pedtools::ped` object, which means the whole
pedsuite is available for building the case: `linearPed`, `cousinPed`,
`avuncularPed`, `halfSibPed`, and profile simulation through `forrel`.

``` r
library(fbnet)
library(pedtools)
library(forrel)

ped <- linearPed(2)                                        # grandparent to grandchild
ped <- setMarkers(ped, locusAttributes = NorwegianFrequencies[1:2])
ped <- profileSim(ped, N = 1, ids = 2, seed = 1)           # genotype the middle generation

pbn  <- from_pedtools(ped, QP = "5")
bn1  <- buildCPTs(buildBN(pbn, QP = "5"))
resQ <- velim.bn(bn1, ordMethod = "min_fill", verbose = FALSE)

names(resQ)
#> [1] "5_D3S1358_p" "5_TH01_p"    "5_D3S1358_m" "5_TH01_m"
```

When `available` is not supplied, the typed individuals are inferred from the
pedigree: anyone carrying at least one non-missing genotype is treated as
evidence.

## The pipeline

| Step | Function | What it does |
|---|---|---|
| 1 | `from_pedtools()` | Reads a `pedtools::ped` into the internal structure: topology, evidence, allele frequencies |
| 2 | `buildBN()` | Builds the network, one node per allele, and marks the query person |
| 3 | `buildCPTs()` | Fills the conditional probability tables. This is the expensive step |
| 4 | `velim.bn()` | Variable elimination, returns the conditioned allele distributions |
| 5 | `reportLR()`, `getQSetRMP()` | Likelihood ratio and random match probability for an observed genotype |

`pruneNodes()`, `stateRemoval()` and `evidencePrunning()` reduce the network
before step 3 and are what keep large pedigrees tractable. `minOrdering()` and
`setOrdering()` control the elimination order used in step 4.

## Legacy entry points

`initBN()`, `initBN.fromPed()`, `initBN.fromVars()`, `convertPedformat()` and
`preparePed()` take pedigrees in the `paramlink` representation. `paramlink` was
archived from CRAN in 2022, so those functions now stop with a message pointing
to `from_pedtools()` unless the package is installed from the CRAN archive.
Existing scripts written against `initBN()` keep working if `paramlink` is
present; new work should use `from_pedtools()`.

## Fragility diagnostics

Per-marker LRs obtained from the network can be passed to `mispitools`
(version 1.4.0 and above) to see how the support is distributed across markers,
and how much of it survives dropping any single one. A case whose LR rests
almost entirely on one marker warrants a different reading from one where the
support is spread out, even when the two report the same total. The package
vignette shows the full workflow.

## References

Chernomoretz, A. et al. (2021). Forensic Science International: Reports.
https://doi.org/10.1016/j.fsir.2020.100132

Darwiche, A. (2009). *Modeling and Reasoning with Bayesian Networks*. Cambridge
University Press.

Vigeland, M. D. (2021). *Pedigree Analysis in R*. Academic Press.

Egeland, T., Kling, D. & Mostad, P. (2015). *Relationship Inference with Familias
and R*. Academic Press.

## Funding

fbnet and GENis are supported by Fundación Sadosky
(https://www.fundacionsadosky.org.ar/) and the Ministerio de Ciencia, Tecnología
e Innovación of Argentina (https://www.argentina.gob.ar/ciencia).
