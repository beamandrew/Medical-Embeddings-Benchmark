library(readr)
library(text2vec)
library(dplyr)
source("./Scripts/utils.R")
options(stringsAsFactors = FALSE)

benchmark_comorbidities <- function(embeddings,k) {
  df <- data.frame(class = character(), concept = character(), dcg = numeric(), associations = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embeddings))){return(NULL)}
  for(file in list.files('./data/benchmarks/comorbidities/')){
    comorbidity <- load_comorbidity(paste0('./data/benchmarks/comorbidities/',file))
    concepts <- intersect(comorbidity$CUI[which(comorbidity$Type=='Concept')],rownames(embeddings))
    missing_concepts <- comorbidity$String[match(setdiff(comorbidity$CUI[which(comorbidity$Type=='Concept')],concepts),comorbidity$CUI)]
    associations <- intersect(comorbidity$CUI[which(comorbidity$Type=='Association')],rownames(embeddings))
    length_df <- dim(df)[1]
    for(i in 1:length(concepts)){
      sim_scores <- get_dist(embeddings,concepts[i])
      df[i+length_df,]<-c(strsplit(file,'.txt'), cuis_to_string(concepts[i]), dcg(sim_scores[1:k], associations), length(associations))
    }
    length_df <- dim(df)[1]
    for(i in 1:length(missing_concepts)){
      df[i+length_df,]<-c(strsplit(file,'.txt'),(missing_concepts[i]),NA,NA)
    }
  }
  return(df)
}

benchmark_semantic_type <- function(embeddings,k){
  df <- data.frame(semantic_type = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embeddings))){return(NULL)}
  for(file in list.files('./data/benchmarks/semantic_type/')){
    length_df <- dim(df)[1]
    semantic_type <- load_semantic_type(paste0('./data/benchmarks/semantic_type/',file))
    cuis <- intersect(semantic_type$CUI,rownames(embeddings))
    if(length(cuis)==0 | k==0){
      df[length_df+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <- 0.0
    for(i in 1:length(cuis)){
      distances <- get_dist(embeddings,cuis[i])[1:k]
      map <- map + length(intersect(cuis[-i],names(distances)))/k
    }
    map <- map / length(cuis)
    df[length_df+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}

benchmark_causitive <- function(embeddings,k){
  df <- data.frame(cause = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embeddings))){return(NULL)}
  for(file in list.files('./data/benchmarks/causative/')){
    length_df <- dim(df)[1]
    causitive <- load_causitive(paste0('./data/benchmarks/causative/',file))
    cuis <- intersect(which(causitive$CUI_Cause %in% rownames(embeddings)),which(causitive$CUI_Result %in% rownames(embeddings)))
    if(length(cuis)==0 | k==0){
      df[length_df+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <- 0.0
    for(i in 1:length(cuis)){
      dist <- get_dist(embeddings,causitive$CUI_Cause[cuis[i]])[1:k]
      if(causitive$CUI_Result[cuis[i]] %in% names(dist)){
        map <- map+1
      }
    }
    map <- map / length(cuis)
    df[length_df+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}

benchmark_ndf_rt <- function(embeddings,k){
  df <- data.frame(ndf_rt = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embeddings))){return(NULL)}
  for(file in list.files('./data/benchmarks/ndf_rt/')){
    length_df <- dim(df)[1]
    ndfrt <- load_ndf_rt(paste0('./data/benchmarks/ndf_rt/',file))
    cuis <- NULL
    for (i in 1:length(ndfrt$Treatment)){
      if((ndfrt$Treatment[i] %in% rownames(embeddings)) & length(intersect(strsplit(ndfrt$Condition[i],','),rownames(embeddings)))>0){
        cuis <- c(cuis,i)
      }
    }
    if(length(cuis)==0 | k==0){
      df[length_df+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <-0.0
    for(i in 1:length(cuis)){
      dist <- get_dist(embeddings, ndfrt$Treatment[cuis[i]])[1:k]
      valid_condition_cui <- intersect(strsplit(ndfrt$Condition[cuis[i]],','),rownames(embeddings))
      map <- map + length(intersect(valid_condition_cui,names(dist)))/length(valid_condition_cui)
    }
    map <- map / length(cuis)
    df[length_df+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}