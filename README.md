
# fbnet

<img src="logo.png" align="right" width="100">

Código con las funciones core de fbnet empaquetado con el criterio de CRAN. Tiene además las distintas funciones comentadas y un manual de referencia. Pasa el test de > R CMD check --as-cran con una nota "New submission". Para usarlo hay que:
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
