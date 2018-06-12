#ensuring github package dependencies available on load
.onLoad <- function(libname, pkgname){
  devtools::install_github("MaayanLab/genesetr")
}