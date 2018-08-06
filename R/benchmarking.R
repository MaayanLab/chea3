benchmarkCDF = function(predict_lib, valid_lib, shuffled_ctl = T,
  background = NULL, method = "FET", integrate = F){
  #check whether valid_lib is nested list
  if(any(sapply(valid_lib, is.list))) error("May only use one validation library.")
  #check whether predict_lib is nested list
  if(any(sapply(predict_lib, is.list))) error("May use one prediction library.
    For multiple prediction libraries use chea3::multibenchmarkCDF()")


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

multibenchmarkCDF = function(predict_libs, valid_lib,
  background = NULL, method = "FET", integrate = T){
  #check whether valid_lib is nested list
  if(any(sapply(valid_lib, is.list))) error("May only use one validation library.")

  #get enrichment results for all prediction libraries
  results = lapply(predict_libs, function(p_lib){
    return(genesetr::pairwiseSetOverlap(lib1 = p_lib, lib2 = valid_lib,
    background = background, method = method))})
  results = lapply(results, function(result) return(chea3::resultsRank(result,2)))
  results = lapply(results, function(result){
    result$TF = unlist(sapply(strsplit(result$set1,"_"),"[",1))
    return(result)
  })
  results = lapply(results, chea3::resultsRank,2)

  if(integrate == T){
    results = chea3::integrateResultsDF(results)
  }

  #shuffle the validation library in each results df
  shuffled_results = lapply(results, function(result){
    return(chea3::shuffleResults(result, shuffle_set = "set2"))})
  shuffled_results = lapply(shuffled_results,
    function(shuffled_result) return(chea3::resultsRank(shuffled_result,"set2")))



  results_cdf = lapply(results,function(result){
    return(ecdf(result[unlist(sapply(strsplit(result$set1,"_"),"[",1)) ==
        unlist(sapply(strsplit(result$set2,"_"),"[",1)),"rank"]))
  })

  shuffled_results_cdf = lapply(shuffled_results, function(shuffled_result){
    return(ecdf(shuffled_result[unlist(sapply(strsplit(shuffled_result$set1,"_"),"[",1)) == unlist(sapply(strsplit(shuffled_result$set2,"_"),"[",1)), "rank"]))
  })

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

