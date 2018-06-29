#' Query ChEA3
#'
#' Compute enrichment of a gene set query over TF regulon libraries
#'
#' @param set_name Name of query. Defaults to "untitled".
#' @param geneset Character array of gene set to query
#' @param n_results Number of results to return for each TF regulon library. Defaults to 10.
#' @param libnames Character array of TF regulon library names to use. Defaults to all libraries available in the package.
#' @param integrate_libs Integrate library results. Often results in a better transcription factor ranking. Defaults to TRUE.
#' @param method Method(s) for comparison of TF regulons and query. Defaults to Fisher's Exact Test
#' @param background The size of the universe of genes under consideration. User may set to a constant. Default is NULL which results in background being set to the size of the union of the unique genes of the query and TF regulon library.
#' @return A list of dataframes with the top TF regulon hits for each TF regulon library tested.
#' @examples
#' data(example_query)
#' results = queryChea3(set_name = "exp1_DEG",example_query)
#'

queryChea3 = function(
  geneset,
  set_name = "untitled",
  n_results = 10,
  libnames = names(chea3::libs),
  integrate_libs = T,
  method = c("FET"),
  background = NULL){

  query = list(geneset)
  names(query) = set_name

  #check library names and ensure in-package availability
  if(!all(libnames %in% names(libs))){
    error(paste("The following libraries are unavailable:",
      paste(setdiff(libnames,names(libs)),collapse = " "),
      "These are the libraries available:",
      paste(names(libs),collapse = " ")))}

  my_libs = chea3::libs[libnames]

  #list of results data frames

  results = lapply(my_libs,function(x){
    df = genesetr::pairwiseSetOverlap(query, x,
      background = background, method = method)
    df$rank = rank(df$FET.p.val,ties.method = "random")/nrow(df)
    df$TF = unlist(sapply(strsplit(df$set2,"_"),"[",1))
    df[,c("a","b","c","d")] = NULL
    return(df)
  })


  if(integrate_libs){
    df_results = plyr::ldply(results,function(sub){
      return(sub)
    })
    int_results = plyr::ddply(df_results,plyr::.(TF),function(sub){
      sub = sub[order(sub$rank),][1,]
      return(sub)
    })
    int_results$rank = rank(int_results$rank,ties.method = "random")

    results[["integrated"]] = int_results

    results = lapply(results,function(sub){
      sub$rank = rank(sub$rank,ties.method = "random")
      return(sub[order(sub$rank),][1:n_results,])
    })
  }

  return(jsonlite::toJSON(results))
}


