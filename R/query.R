#' Query ChEA3
#'
#' Compute enrichment of a gene set query over TF regulon libraries
#'
#' @param set_name Name of query. Defaults to "untitled".
#' @param geneset Character array of gene set to query
#' @param n_results Number of results to return for each TF regulon library. Defaults to 10. Input "all"
#' to return results for all transcription factors.
#' @param libnames Character array of TF regulon library names to use. Defaults to all libraries available in the package.
#' @param integrate_libs Integrate library results. Often results in a better transcription factor ranking. Defaults to TRUE.
#' @param method Method(s) for comparison of TF regulons and query. Defaults to Fisher's Exact Test
#' @param background The size of the universe of genes under consideration. User may set to a constant. Default is NULL which results in background being set to the size of the union of the unique genes of the query and TF regulon library.
#' @return A list of dataframes with the top TF regulon hits for each TF regulon library tested.
#' @examples
#' data(example_query)
#' results = queryChea3(set_name = "exp1_DEG",example_query)
#'

queryChea = function(
  geneset,
  set_name = "untitled",
  n_results = 10,
  libnames = names(chea3::libs),
  integrate_libs = T,
  integrate_libnames = setdiff(names(chea3::libs),"BioGRID"),
  method = c("FET"),
  background = NULL){

  query <- list(geneset)
  names(query) <- set_name

  #check library names and ensure in-package availability
  if(!all(libnames %in% names(chea3::libs))){
    error(paste("The following libraries are unavailable:",
      paste(setdiff(libnames,names(libs)),collapse = " "),
      "These are the libraries available:",
      paste(names(libs),collapse = " ")))}

  my_libs <- chea3::libs[libnames]

  #list of results data frames

  results <- lapply(my_libs,function(x){
    x <- genesetr::removeLibWeights(x)
    df <- genesetr::pairwiseSetOverlap(x, query,
      background = background, method = method)
    df$rank <- rank(df$FET.p.val,ties.method = "random")/nrow(df)
    df$TF <- unlist(sapply(strsplit(df$set1,"_"),"[",1))
    df$FET.p.val <- signif(df$FET.p.val,3)
    df[,c("a","b","c","d")] <- NULL
    df <- genesetr::dfcols.tochar(df)
    return(df)
  })

  if(integrate_libs){
    results <- chea3::integrateResultsDF(results,integrate_libnames)}

  results <- lapply(results, function(result){
    result$rank <- rank(result$rank, ties.method = "random")
    result <- result[order(result$rank),]
    if(n_results == "all"){
      return(result)
      }
    else{
        return(result[1:n_results,])
    }
  })
  return(results)
}



#' Web query ChEA3
#'
#' Serves as a wrapper function for queryChea()
#'
#' @param set_name Name of query. Defaults to "untitled".
#' @param geneset Character array of gene set to query
#' @param n_results Number of results to return for each TF regulon library. Defaults to 10. Input "all"
#' to return results for all transcription factors.
#' @param libnames Character array of TF regulon library names to use. Defaults to all libraries available in the package.
#' @param integrate_libs Integrate library results. Often results in a better transcription factor ranking. Defaults to TRUE.
#' @param method Method(s) for comparison of TF regulons and query. Defaults to Fisher's Exact Test
#' @param background The size of the universe of genes under consideration. User may set to a constant. Default is NULL which results in background being set to the size of the union of the unique genes of the query and TF regulon library.
#' @return A json of datasets with the top TF regulon hits for each TF regulon library tested.
#' @examples
#' data(example_query)
#' results = queryChea3(set_name = "exp1_DEG",example_query)
#'
queryCheaWeb = function(
  geneset,
  set_name = "untitled",
  n_results = "all",
  libnames = names(chea3::libs),
  integrate_libs = T,
  method = c("FET"),
  background = 30000){

  results <- queryChea(geneset = geneset,
    set_name = set_name,
    n_results = n_results,
    libnames = libnames,
    integrate_libs = integrate_libs,
    method = method,
    background = background)

  # lib_descripts <- chea3::getLibDescriptShort()
  # lib_descripts <- lib_descripts[order(unlist(lib_descripts))]
  # lib_descripts <- rlist::list.prepend(lib_descripts,Integrated = "integrated results from all resources")
  #
  # results <- results[match(names(lib_descripts),names(results))]
  # idx <- match(names(results),names(lib_descripts))
  #
  # names(results) <- paste(names(results), ": ",
  #   unlist(lib_descripts[idx]),sep = "")

  json_results =  gsub("FET.p.val","FET p-value",jsonlite::toJSON(results))

  return(jsonlite::toJSON(results))
}

integrateResultsDF = function(results, integrate_libnames){
  # if(!sum(sapply(list,is.data.frame))>1) error("Must have >1 results dataframe in order to integrate results.")

  df_results <- plyr::ldply(results,function(sub){
    return(sub)
  })
  df_results <- df_results[df_results[,".id"] %in% integrate_libnames,]
  int_results <- plyr::ddply(df_results,plyr::.(set2,TF),function(sub){
    sub = sub[order(sub$rank),][1,]
    return(sub)
  })
  int_results$prev_rank <- int_results$rank
  int_results$rank <- rank(int_results$rank, ties.method = "random")/nrow(int_results)
  int_results$.id <- as.character(int_results$.id)

  results[["Integrated"]] <- int_results
  return(results)
} #end function integrateResultsDF()
