<img src="logo.png" align="right" width="120">

# fbnet
fbnet is an open source software package written in R statistical languaje. It relies on a Bayesian Networks framework and it is particularly well suited to efficiently perform large-size queries against databases of missing individuals.

It could interact with the main functionallities of other packages for pedigree analysis. In particular, fbnet can interpret the pedigree formats from 'Familias' software. In addition 'pedtools' format, a software for creating and manipulating pedigrees and markers, is supported. 

fbnet allows computing LRs and obtaining genotype probability distributions for query individuals, based on the pedigree data. 

fbnet implements the complete GENis functionality, a recently published open-source multi-tier information system developed to run forensic DNA databases to perform kinship analysis based on DNA profiles. For more information about GENis see [Chernomoretz et al.](https://www.sciencedirect.com/science/article/pii/S2665910720300815?via%3Dihub)
.

Stable version could be installed from CRAN:

      > install.packages("fbnet") 
      > library(fbnet)


Version under development could be installed as follows: 

      > library(devtools)
      > install_github("MarsicoFL/fbnet")
      > library(fbnet)
      
 Example data could be parsed as follows:
 
      > pbn  <- initBN(toyped)
      > bnet <- buildBN(pbn,QP=3)
      > bn1  <- buildCPTs(bnet)
      > resQ <- velim.bn(bn1,ordMethod="min_fill",verbose=FALSE)
