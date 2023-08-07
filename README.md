<img src="genis.png" align="left" width="80">



# fbnet

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/fbnet)](https://CRAN.R-project.org/package=fbnet)
[![](https://cranlogs.r-pkg.org/badges/grand-total/fbnet?color=blue)](https://cran.r-project.org/package=fbnet)

<!-- badges: end -->

## IMPORTANT NOTE
This package is only maintained for legacy purposes. New releases of the fbnet package and under-development version are available on https://github.com/chernolabs/fbnet-cran.


fbnet is an open source software package written in R statistical languaje. It implements the complete *GENis* functionality, a recently published open-source multi-tier information system developed to run forensic DNA databases to perform kinship analysis based on DNA profiles. GENis is freely available on github: https://github.com/fundacion-sadosky/genis


It relies on a Bayesian Networks framework and it is particularly well suited
to efficiently perform large-size queries against databases of missing individuals.
It could interact with the main functionallities of other packages for pedigree analysis. 
In particular, fbnet can interpret the pedigree formats from 'Familias' software (1). In addition 'pedtools', a software for creating and manipulating pedigrees and markers, is supported (2). fbnet allows computing LRs
and obtaining genotype probability distributions for a query individual, based on 
the pedigree data. For installing fbnet please use the following command:

``` r
install.packages("fbnet")
library(fbnet)
```

fbnet and GENis projects are supported by Fundación Sadosky https://www.fundacionsadosky.org.ar/ and Ministerio de Ciencia, Tecnología e Innovación of Argentina https://www.argentina.gob.ar/ciencia


1- Vigeland, M. D. (2021). Pedigree analysis in R. Academic Press. <br /> 
2- Egeland, T., Kling, D., & Mostad, P. (2015). Relationship inference with familias and R: statistical methods in forensic genetics. Academic Press.
