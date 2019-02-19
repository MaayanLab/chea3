#package_data
libs = chea3::libs
libs["BioGRID"] = NULL

libs[["Perturbations"]] = genesetr::loadGMT("/Users/alexandrakeenan/Desktop/USE THIS FOR SINGLE-TF PERTURBATIONS/singleTFperturbs.gmt")


libs = lapply(libs, genesetr::removeDupes)
libs = lapply(libs, genesetr::removeEmptySets)

devtools::use_data(libs, internal = T, overwrite = T)


Perturbations = libs[["Perturbations"]]
devtools::use_data(Perturbations, overwrite = T)

