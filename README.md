<img src="logo.png" align="right" width="120">

# fbnet
fbnet is an open source software package written in R statistical languaje.
It relies on a Bayesian Networks framework and it is particularly well suited
to efficiently perform large-size queries against databases of missing individuals.
It could interact with the main functionallities of other packages for pedigree analysis. 
In particular, fbnet can interpret the pedigree formats from 'Familias' software. In addition 'pedtools', a software for creating and manipulating pedigrees and markers, is supported. fbnet allows computing LRs
and obtaining genotype probability distributions for query individual, based on 
the pedigree data. fbnet implements the complete GENis functionality, a recently published open-source 
multi-tier information system developed to run forensic DNA databases  
to perform kinship analysis based on DNA profiles.
Para usarlo hay que:



1) descargar el código (click en code y luego dowload).
2) Descomprimirlo
3) Abrir el archivo fbnet.Rproj con R o Rstudio.
4) Instalar y cargar devtools.
      > install.packages("devtools") 
      
      > library(devtools)
5) Ejectutar algunas funciones de devtools para chequear instalación.
      > load_all()

      > document()

      > install()
6) Cargar fbnet y testear funciones
      > library(fbnet)

      > ?fbnet()
