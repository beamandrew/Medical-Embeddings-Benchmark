library(readr)
library(text2vec)
library(dplyr)
source("./Scripts/utils.R")
options(stringsAsFactors = FALSE)

benchmark_comorbidities <- function(embedding,k, ref_cuis=NULL, return_max = FALSE) {
  df <- data.frame(class = character(), concept = character(), dcg = numeric(), associations = numeric(), stringsAsFactors = FALSE)
  
  for(file in list.files('./data/benchmarks/comorbidities/')){
    comorbidity <- load_comorbidity(paste0('./data/benchmarks/comorbidities/',file))
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    if(k>min(length(rownames(embedding)),length(ref_cuis))){return(NULL)}
    concepts <- intersect(comorbidity$CUI[which(comorbidity$Type=='Concept')],ref_cuis)
    strings <- comorbidity$String[comorbidity$CUI %in% concepts]
    associations <- intersect(comorbidity$CUI[which(comorbidity$Type=='Association')],ref_cuis)
    max_score <- 0 
    max_string <- ''
    for(i in 1:length(concepts)){
      sim_scores <- get_dist(embedding,concepts[i])[1:k]
      score <- dcg(sim_scores, associations)
      if(!return_max){
        df[dim(df)[1]+1,]<-c(strsplit(file,'.txt'), strings[i], score, length(associations))
      }
      if(return_max){
        if(score>max_score){
          max_score <- score
          max_string <- strings[i]
        }
      }
    }
    if(return_max){
      df[dim(df)[1]+1,]<-c(strsplit(file,'.txt'),'max',max_score,length(associations))
    }
  }
  return(df)
}

benchmark_semantic_type <- function(embedding,k,ref_cuis=NULL){
  df <- data.frame(semantic_type = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embedding))){return(NULL)}
  for(file in list.files('./data/benchmarks/semantic_type/')){
    length_df <- dim(df)[1]
    semantic_type <- load_semantic_type(paste0('./data/benchmarks/semantic_type/',file))
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    cuis <- intersect(semantic_type$CUI,ref_cuis)
    if(length(cuis)==0 | k==0){
      df[length_df+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <- 0.0
    for(i in 1:length(cuis)){
      distances <- get_dist(embedding,cuis[i])[1:k]
      map <- map + length(intersect(cuis[-i],names(distances)))/k
    }
    map <- map / length(cuis)
    df[length_df+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}

benchmark_causitive <- function(embedding,k,ref_cuis=NULL){
  df <- data.frame(cause = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embedding))){return(NULL)}
  for(file in list.files('./data/benchmarks/causative/')){
    length_df <- dim(df)[1]
    causitive <- load_causitive(paste0('./data/benchmarks/causative/',file))
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    cuis <- intersect(which(causitive$CUI_Cause %in% ref_cuis),which(causitive$CUI_Result %in% ref_cuis))
    if(length(cuis)==0 | k==0){
      df[length_df+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <- 0.0
    for(i in 1:length(cuis)){
      dist <- get_dist(embedding,causitive$CUI_Cause[cuis[i]])[1:k]
      if(causitive$CUI_Result[cuis[i]] %in% names(dist)){
        map <- map+1
      }
    }
    map <- map / length(cuis)
    df[length_df+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}

benchmark_ndf_rt <- function(embedding,k,ref_cuis=NULL){
  df <- data.frame(ndf_rt = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embedding))){return(NULL)}
  for(file in list.files('./data/benchmarks/ndf_rt/')){
    length_df <- dim(df)[1]
    ndfrt <- load_ndf_rt(paste0('./data/benchmarks/ndf_rt/',file))
    cuis <- NULL
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    for (i in 1:length(ndfrt$Treatment)){
      if((ndfrt$Treatment[i] %in% ref_cuis) & length(intersect(strsplit(ndfrt$Condition[i],','),ref_cuis))>0){
        cuis <- c(cuis,i)
      }
    }
    if(length(cuis)==0 | k==0){
      df[length_df+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <-0.0
    for(i in 1:length(cuis)){
      dist <- get_dist(embedding, ndfrt$Treatment[cuis[i]])[1:k]
      valid_condition_cui <- intersect(strsplit(ndfrt$Condition[cuis[i]],','),ref_cuis)
      map <- map + length(intersect(valid_condition_cui,names(dist)))/length(valid_condition_cui)
    }
    map <- map / length(cuis)
    df[length_df+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}