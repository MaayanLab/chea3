libs2network = function(libnames = chea3::getLibNames()){
  edges = data.frame()
  for(i in 1:length(libnames)){
    lib = getLib(libnames[i])
    lib = genesetr::removeDupes(lib)
    lib = genesetr::toLongDF(lib)
    lib$source = unlist(sapply(strsplit(lib$set_name,"_"),"[",1))
    edges = rbind(edges,data.frame(source = lib$source, target = lib$gene, lib = libnames[i]))
  }
  return(edges)
}