library(dplyr)
library(readr)
library(text2vec)
library(magrittr)
id_file <- read_delim("./data/2a_concept_ID_to_string.txt",delim="\t",col_names = FALSE,quote="")
colnames(id_file) <- c("Concept_ID","String")
cui_file <- read_delim("./data/2b_concept_ID_to_CUI.txt",delim="\t",col_names = FALSE,quote="")
colnames(cui_file) <- c("Concept_ID","CUI")

info <- id_file %>% inner_join(cui_file)

#' Compute the mean average precision at k
#'
#' This function computes the mean average precision at k
#' of two lists of sequences.
#'
#' @param k max length of predicted sequence
#' @param actual list of ground truth sets (vectors)
#' @param predicted list of predicted sequences (vectors)
#' @export
#' Adapted from Kaggle: https://github.com/benhamner/Metrics/blob/master/R/R/metrics.r#L181
mapk <- function (k, actual, predicted)
{
  if( length(actual)==0 || length(predicted)==0 ) 
  {
    return(0.0)
  }
  
  scores <- rep(0, length(actual))
  for (i in 1:length(scores))
  {
    scores[i] <- apk(k, actual[[i]], predicted[[i]])
  }
  score <- mean(scores)
  score
}


#' Compute the average precision at k
#'
#' This function computes the average precision at k
#' between two sequences
#'
#' @param k max length of predicted sequence
#' @param actual ground truth set (vector)
#' @param predicted predicted sequence (vector)
#' @export
apk <- function(k, actual, predicted)
{
  score <- 0.0
  cnt <- 0.0
  for (i in 1:min(k,length(predicted)))
  {
    if (predicted[i] %in% actual && !(predicted[i] %in% predicted[0:(i-1)]))
    {
      cnt <- cnt + 1
      score <- score + cnt/i 
    }
  }
  score <- score / min(length(actual), k)
  score
}

get_dist <- function(word_vectors,query,sort_result=TRUE) {
  word_vectors <- as.matrix(word_vectors)
  word_vectors_norm <- sqrt(rowSums(word_vectors ^ 2))
  query_vec <- word_vectors[query,,drop=FALSE]
  cos_dist <- text2vec:::cosine(query_vec, 
                                word_vectors, 
                                word_vectors_norm)
  
  cos_dist <- cos_dist[1,]
  if(sort_result) {
    cos_dist <- sort(cos_dist,decreasing = TRUE)
  }
  return(cos_dist)
}

cuis_to_string <- function(cuis) {
  strings <- info$String[match(cuis,info$CUI)]
  return(strings)
}

load_embeddings <- function(filename,convert_to_cui=TRUE,header=F,skip=1) {
  embeddings <- read.delim(filename,sep=" ",skip=skip,header=header)
  rownames(embeddings) <- embeddings[,1]
  embeddings <- embeddings[,-1]
  if(convert_to_cui) {
    cuis <- info$CUI[match(rownames(embeddings),info$Concept_ID)]
    rownames(embeddings) <- cuis
  }
  return(embeddings)
}


dcg <- function(vector, true_list){
  score <- 0 
  cuis <- names(vector)
  relevant_cuis <- which(cuis %in% true_list)
  if (length(relevant_cuis)==0){return(0)}
  for(i in 1:length(relevant_cuis)){
    score <- score + (2^vector[relevant_cuis[i]]-1)/log2(relevant_cuis[i])
  }
  return(score)
}

load_comorbidity <- function(filename){
  commorbidity <- read.delim(filename)
  return(commorbidity)
}

load_semnatic_type <- function(filename,header=F,skip=1) {
  semantic <- read.delim(filename)
  return(semantic)
}




