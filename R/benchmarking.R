benchmarkCDF = function(predict_lib, valid_lib, shuffled_ctl = T,
  background = NULL, method = "FET"){
  results = genesetr::pairwiseSetOverlap(lib1 = predict_lib, lib2 = valid_lib,
    background = background, method = method)
  #shuffle the validation library
  shuffled_results = chea3::shuffleResults(results, shuffle_set = 2)

  results = chea3::resultsRank(results,2)
  shuffled_results = chea3::resultsRank(shuffled_results,2)

  results_cdf = ecdf(results[unlist(sapply(strsplit(results$set1,"_"),"[",1)) ==
      unlist(sapply(strsplit(results$set2,"_"),"[",1)),"rank"])

  shuffled_results_cdf = ecdf(shuffled_results[unlist(sapply(strsplit(shuffled_results$set1,"_"),"[",1)) == unlist(sapply(strsplit(shuffled_results$set2,"_"),"[",1)), "rank"])

  return(list(results_cdf = results_cdf, shuffled_results_cdf = shuffled_results_cdf))
}

shuffleResults = function(results,shuffle_set = 2){
  names = unique(results[,shuffle_set])
  shuffled_names = sample(names,replace = F)
  shuffled_results = results
  for(i in 1:length(names)){
    idx = results[,shuffle_set] == names[i]
    shuffled_results[idx,shuffle_set] = shuffled_names[i]
  }
  return(shuffled_results)
}

resultsRank = function(results, valid_set){
  return(plyr::ddply(results,valid_set,function(sub){
    sub$rank = rank(sub$FET.p.val,ties.method = "random")/nrow(sub)
    return(sub)
  }))
}

