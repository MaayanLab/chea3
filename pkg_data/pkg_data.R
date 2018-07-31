#package_data
libs = list()
libs[["ReMap"]] = genesetr::loadGMT("/Users/alexandrakeenan/Projects/TF_libs/gmts/all_remap/0.05_remap_all.tsv")
libs[["ARCHS4"]] = genesetr::loadGMT("/Users/alexandrakeenan/Projects/TF_libs/gmts/coexpression/ARCHS4_TFs_pearson_n_50000_human.gmt")
libs[["ChEA"]] = genesetr::loadGMT("/Users/alexandrakeenan/Projects/TF_libs/gmts/enrichr_dwnld_derived/chipseq/HGNC_mapped_enrichr_ChEA_2016.txt")
libs[["GTEx"]] = genesetr::loadGMT("/Users/alexandrakeenan/Projects/TF_libs/gmts/coexpression/GTEx_TFs_pearson.gmt")
libs[["Enrichr"]] = genesetr::loadGMT("/Volumes/External/CHEA3/Enrichr_TF_cooccur/Enrichr_tf_cooccur_probability.gmt")

libs = lapply(libs, genesetr::removeDupes)
libs = lapply(libs, genesetr::removeEmptySets)

devtools::use_data(libs, internal = T, overwrite = T)

ARCHS4 = libs[["ARCHS4"]]
devtools::use_data(ARCHS4, overwrite = T)

ChEA = libs[["ChEA"]]
devtools::use_data(ChEA, overwrite = T)

ReMap = libs[["ReMap"]]
devtools::use_data(ReMap, overwrite = T)

Enrichr = libs[["Enrichr"]]
devtools::use_data(Enrichr, overwrite = T)

GTEx = libs[["GTEx"]]
devtools::use_data(GTEx, overwrite = T)