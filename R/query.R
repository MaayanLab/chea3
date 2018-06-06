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

queryChea3 = function(set_name = "untitled",
  geneset,
  n_results = 10,
  libnames = names(libs),
  integrate_libs = T,
  method = c("FET"),
  background = NULL){

  query = list(geneset)
  names(query) = set_name


  if(!all(libnames %in% names(libs))){
    error(paste("The following libraries are unavailable:",
      paste(setdiff(libnames,names(libs)),collapse = " "),
      "These are the libraries available:",
      paste(names(libs),collapse = " ")))}

  temp_libs = libs[libnames]
  results = lapply(temp_libs,function(x){
    df = genesetr::pairwiseSetOverlap(query, x,
      background = background, method = method)
    df = df[order(df$FET.p.val),][1:n_results,]
    return(df)
  })
  return(results)
}



