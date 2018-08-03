library(magrittr)
stats = chea3::libs %>% lapply(genesetr::getLibInfo) %>% plyr::ldply(data.frame,stringsAsFactors = F)


stats = plyr::ddply(stats, plyr::.(.id), function(stat){
  set_names = names(chea3::libs[[stat$.id]])
  stat$n_unique_tfs = set_names %>% strsplit("_") %>% lapply(head,1) %>%
    unique %>% length
  return(stat)
})

write.table(stats,'/users/alexandrakeenan/Desktop/Vivian_lib_stats.tsv',sep = "\t",col.names = T, row.names = F, quote = F)
