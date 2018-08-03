lib_meta = list()

lib_meta[["ARCHS4"]] = "Co-expression; TF-gene co-expression from uniformly reprocessed GEO RNA-seq experiments"
lib_meta[["GTEx"]] = "Co-expression; TF-gene co-epxression from GTEx RNA-seq data"
lib_meta[["ChEA"]] = "ChIP-seq; Manually curated TF regulons from publications with ChIP-seq experiments"
lib_meta[["ReMap"]] = "ChIP-seq; Uniformly re-processed ChIP-seq experiments from GEO and ENCODE"
lib_meta[["ENCODE"]] = "ChIP-seq; TF targets from the ENCODE consortium"
lib_meta[["Enrichr"]] = "Co-occurrence; TF-gene co-occurrence in crowd-submitted gene sets"
lib_meta[["Perturbations"]] = "Gene signatures; Differentially expressed genes on genetic perturbation of a single TF (e.g. KO, KD, OE)"
lib_meta[["BioPLEX"]] = "Protein-protein interactions; TF protein interactions curated from published studies"

devtools::use_data(lib_meta, overwrite = T)
