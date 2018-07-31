libs2network = function(libnames = chea3::getLibNames()){
  edges = data.frame()
  for(i in 1:length(libnames)){
    lib = chea3::getLib(libnames[i])
    lib = genesetr::removeDupes(lib)
    lib = genesetr::toLongDF(lib)
    lib$source = unlist(sapply(strsplit(lib$set_name,"_"),"[",1))
    edges = rbind(edges,data.frame(source = lib$source, target = lib$gene, lib = libnames[i], weight = lib$weight, stringsAsFactors = F))
  }
  return(edges)
}

toTFNet = function(network){
  return(network[network$source %in% chea3::tfs &
      network$target %in% chea3::tfs, ])
}