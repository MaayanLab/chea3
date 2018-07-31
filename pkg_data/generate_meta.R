lib_meta = list()
lib_meta[["ReMap"]] = "ChIP-seq; Uniformly re-processed ChIP-seq experiments from ENCODE and GEO"
lib_meta[["ARCHS4"]] = "Co-expression; TF-gene co-expression from uniformly reprocessed GEO RNA-seq experiments"
lib_meta[["ChEA"]] = "ChIP-seq; Manually curated TF regulons from publications with ChIP-seq experiments"
lib_meta[["GTEx"]] = "Co-expression; TF-gene co-epxression from GTEx RNA-seq data"
lib_meta[["Enrichr"]] = "Co-occurrence; TF-gene co-occurrence in crowd-submitted gene sets"

devtools::use_data(lib_meta)
