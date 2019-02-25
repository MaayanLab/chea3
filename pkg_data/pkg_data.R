#package_data
rm(list = ls())

#prepare TFs data
lamb_tfs = read.table("/volumes/external/Projects/Common/Lambert2018_TFs/formatted_human_tfs.txt",stringsAsFactors=F, quote="", comment.char="", sep="\t")$V2
mapped_lamb_tfs = genesetr::HGNCapproved(lamb_tfs, untranslatable.na = T)
mapped_lamb_tfs = unique(mapped_lamb_tfs[!is.na(mapped_lamb_tfs)])
tfs = mapped_lamb_tfs
devtools::use_data(tfs, overwrite = T)


libs = list()

libs[["Perturbations"]] = genesetr::loadGMT("/Volumes/External/Projects/TF_libs/gmts/single_TF_pert/RNAseq_Micro_TF_perts_all.gmt")
libs[["ChEA"]] = genesetr::loadGMT("/Volumes/External/Projects/TF_libs/gmts/chea_encode_chipseq_gmts/ChEA2015_mapped.gmt")
libs[["ENCODE"]] = genesetr::loadGMT("/Volumes/External/Projects/TF_libs/gmts/chea_encode_chipseq_gmts/ENCODE2015_mapped.gmt")
libs[["GTEx"]] = genesetr::loadGMT("/Volumes/External/Projects/TF_libs/gmts/coexpression_v2/GTEx_TFs_pearson.gmt")
libs[["ARCHS4"]] = genesetr::loadGMT("/Volumes/External/Projects/TF_libs/gmts/coexpression_v2/ARCHS4_TFs_pearson_n_50000_human.gmt")
libs[["Enrichr"]] = genesetr::loadGMT("/Volumes/External/Projects/TF_libs/gmts/enricher_cooccurrence/ enrichr_coocurrence.gmt")
libs[["ReMap"]] = genesetr::loadGMT("/Volumes/External/Projects/TF_libs/gmts/all_remap/0.05_remap_all.tsv")

libs = lapply(libs, genesetr::removeLibWeights)
libs = lapply(libs, genesetr::removeDupes)
libs = lapply(libs, genesetr::removeEmptySets)

devtools::use_data(libs, internal = T, overwrite = T)


Perturbations = libs[["Perturbations"]]
devtools::use_data(Perturbations, overwrite = T)

ENCODE = libs[["ENCODE"]]
devtools::use_data(ENCODE, overwrite = T)

ChEA = libs[["ChEA"]]
devtools::use_data(ChEA, overwrite = T)

ReMap = libs[["ReMap"]]
devtools::use_data(ReMap, overwrite = T)

ARCHS4 = libs[["ARCHS4"]]
devtools::use_data(ARCHS4, overwrite = T)

Enrichr = libs[["Enrichr"]]
devtools::use_data(Enrichr, overwrite = T)

GTEx = libs[["GTEx"]]
devtools::use_data(GTEx, overwrite = T)

#write to folder for use in web service
genesetr::writeGMT(GTEx,"/users/alexandrakeenan/eclipse-workspace/chea3/WebContent/WEB-INF/tflibs/GTEx--Coexpression_fromRpackage.txt")
genesetr::writeGMT(ARCHS4,"/users/alexandrakeenan/eclipse-workspace/chea3/WebContent/WEB-INF/tflibs/ARCHS4--Coexpression_fromRpackage.txt")
genesetr::writeGMT(ChEA,"/users/alexandrakeenan/eclipse-workspace/chea3/WebContent/WEB-INF/tflibs/Literature--ChIP-seq_fromRpackage.txt")
genesetr::writeGMT(ENCODE,"/users/alexandrakeenan/eclipse-workspace/chea3/WebContent/WEB-INF/tflibs/ENCODE--ChIP-seq_fromRpackage.txt")
genesetr::writeGMT(ReMap,"/users/alexandrakeenan/eclipse-workspace/chea3/WebContent/WEB-INF/tflibs/ReMap--ChIP-seq_fromRpackage.txt")
genesetr::writeGMT(Enrichr,"/users/alexandrakeenan/eclipse-workspace/chea3/WebContent/WEB-INF/tflibs/Enrichr--Queries_fromRpackage.txt")

