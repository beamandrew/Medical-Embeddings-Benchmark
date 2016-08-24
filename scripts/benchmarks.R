library(readr)
library(text2vec)
library(dplyr)
source("./Scripts/utils.R")
options(stringsAsFactors = FALSE)

#Benchmarking comorbidities takes in text files with Concepts and Associations and computes either the DCG or MAP for each concept
#Return max returns the max DCG over all concepts for one file 
benchmark_comorbidities <- function(embedding,k, ref_cuis=NULL, return_max = FALSE, metric='DCG') {
  #Create the data frame that we will append all scores to 
  df <- data.frame(class = character(), concept = character(), score = numeric(), associations = numeric(), stringsAsFactors = FALSE)
  #Looping over all the comorbidities
  for(file in list.files('./data/benchmarks/comorbidities/')){
    #Load the file
    comorbidity <- load_comorbidity(paste0('./data/benchmarks/comorbidities/',file))
    #If there are no reference CUIs, we need to take the rownames of the embedding as our reference (ie there is no intersection)
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    if(k>min(length(rownames(embedding)),length(ref_cuis))){return(NULL)}
    #Get only the concepts and associations that are in the embedding 
    concepts <- intersect(comorbidity$CUI[which(comorbidity$Type=='Concept')],ref_cuis)
    strings <- comorbidity$String[comorbidity$CUI %in% concepts]
    associations <- intersect(comorbidity$CUI[which(comorbidity$Type=='Association')],ref_cuis)
    max_score <- 0 
    max_string <- ''
    for(i in 1:length(concepts)){
      #Query the embedding and get the top K vectors back
      sim_scores <- get_dist(embedding,concepts[i])[1:k]
      #Compute the DCG 
      if(metric=='DCG'){
        score <- dcg(sim_scores, associations)
      }
      #Compute the Precision 
      if(metric=='AP'){
        score <- length(intersect(names(sim_scores),associations))/k
      }
      #If we are not return the max, then we need to append every score to our returned data frame 
      if(!return_max){
        df[dim(df)[1]+1,]<-c(strsplit(file,'.txt'), strings[i], score, length(associations))
      }
      #Otherwise, we need to keep track of the max 
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

#' Compute the mean average precision at k for 4 semantic type benchmarks derived from the UMLS: cellular_molecular_dysfunction, genetic_function, mental_or_behavioral,
#' and neoplastic_process 
#'
#' 
#' @param embedding Matrix containing the embeddings
#' @param k Cutoff to use for in MAPk
#' @param ref_cuis CUIs to consider. Useful if you are comparing embeddings from different sources that have some non-overlapping concepts.
#' @param verbose Should information be printed to the screen during the benchmark?
#' @export
benchmark_semantic_type <- function(embedding,k,ref_cuis=NULL,verbose=TRUE){
  #Generating the data frame we will return 
  df <- data.frame(semantic_type = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embedding))){return(NULL)}
  #Looping over all the semantic types
  for(file in list.files('./data/benchmarks/semantic_type/')){
    #Load the file
    if(verbose){
      print(paste0("Now benchmarking semantic type: ",strsplit(file,split="\\.")[[1]][1])) 
    }
    semantic_type <- load_semantic_type(paste0('./data/benchmarks/semantic_type/',file))
    #If there are no reference CUIs, there is no intersection taken, so we just take the rownames of the embedding
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    #Need to get only the CUIs that are in the embedding 
    cuis <- intersect(semantic_type$CUI,ref_cuis)
    if(length(cuis)==0 | k==0){
      df[dim(df)[1]+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <- 0.0
    if(verbose) {
      pb <- txtProgressBar(min=1,max=length(cuis),style=3)
    }
    for(i in 1:length(cuis)){
      if(verbose) {
        setTxtProgressBar(pb,i)
      }
      num_hits = 0
      #Querying the embedding for the top K vectors for a CUI
      distances <- get_dist(embedding,cuis[i],sort_result = TRUE)
      ap <- apk(k,cuis[-i],names(distances)[-1])
      map <- map + ap
    }
    # Diving by the number of terms gives us MAP 
    map <- map / length(cuis)
    print(paste0("MAP: ",map))
    df[dim(df)[1]+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}

#Benchmark causitive takes in a file of cause:result pairs and computes the 
#accuracy of having the result in the top k vectors returned when 
#an embedding is queried with the corresponding cause 
benchmark_causitive <- function(embedding,k,ref_cuis=NULL,verbose = TRUE){
  #Generating the data frame we will return 
  df <- data.frame(cause = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embedding))){return(NULL)}
  #Looping over all the causitive files 
  for(file in list.files('./data/benchmarks/causative/')){
    if(verbose){
      print(paste0("Now benchmarking semantic type: ",strsplit(file,split="\\.")[[1]][1])) 
    }
    causitive <- load_causitive(paste0('./data/benchmarks/causative/',file))
    #If there are no reference CUIs we simply take all the CUIs in the embedding 
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    #Need to get only the cause:result pairs where both cause AND result are in the embedding 
    rows_in_both <- intersect(which(causitive$CUI_Cause %in% ref_cuis),which(causitive$CUI_Result %in% ref_cuis))
    cuis <- causitive[rows_in_both,]
    if(nrow(cuis)==0 | k==0){
      df[dim(df)[1]+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <- 0.0
    if(verbose) {
      pb <- txtProgressBar(min=1,max=nrow(cuis),style=3)
    }
    for(i in 1:nrow(cuis)){
      if(verbose) {
        setTxtProgressBar(pb,i)
      }
      #Query the embedding 
      dist <- get_dist(embedding,cuis$CUI_Cause[i],sort_result = TRUE)
      map <- map + apk(k,cuis$CUI_Result[i],names(distances)[-1])
      
    }
    #Dividing by the number of terms gives MAP 
    map <- map / length(cuis)
    print(paste0("MAP: ",map))
    #Writing the score to the data frame 
    df[dim(df)[1]+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}

#Benchmark NDF RT benchmarks treatment/prevention files as determined by the National 
#Drug File - Reference Terminology of the National Library of Medicine
#It is similar to benchmark causitive but now instead of 1
benchmark_ndf_rt <- function(embedding,k,ref_cuis=NULL){
  #Generate the data frame we will return 
  df <- data.frame(ndf_rt = character(), map = numeric(), stringsAsFactors = FALSE)
  if(k>length(rownames(embedding))){return(NULL)}
  #Loop over all the NDF RT files
  for(file in list.files('./data/benchmarks/ndf_rt/')){
    ndfrt <- load_ndf_rt(paste0('./data/benchmarks/ndf_rt/',file))
    #Again, if no reference CUIs, need to use all the CUIs in the embedding
    if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
    #Need to get a list of all treatment cuis with at least one cui in the condition column. Initializing this list as null 
    cuis <- NULL
    for (i in 1:length(ndfrt$Treatment)){
      #If the treatment cui is in the embedding and has at least one condition cui, then we include this treatment cui in the list of valid cuis
      if((ndfrt$Treatment[i] %in% ref_cuis) & length(intersect(strsplit(ndfrt$Condition[i],','),ref_cuis))>0){
        #Appending the valid cui to the list 
        cuis <- c(cuis,i)
      }
    }
    if(length(cuis)==0 | k==0){
      df[dim(df)[1]+1,] <- c(strsplit(file,'.txt'),0)
      next
    }
    map <-0.0
    #Loop over all the valid treatment cuis 
    for(i in 1:length(cuis)){
      #Query the embedding with a valid treatment cui 
      dist <- get_dist(embedding, ndfrt$Treatment[cuis[i]])[1:k]
      #Get all the valid condition cuis for the embedding, guaranteed to have length of at least one since this is a valid treatment CUI 
      valid_condition_cui <- intersect(strsplit(ndfrt$Condition[cuis[i]],','),ref_cuis)
      #Adding the average precision for each treatment CUI
      map <- map + length(intersect(valid_condition_cui,names(dist)))/length(valid_condition_cui)
    }
    #Computing MAP 
    map <- map / length(cuis)
    #Writing the result to the data frame
    df[dim(df)[1]+1,] <- c(strsplit(file,'.txt'),map)
  }
  return(df)
}

#Benchmark similarity computes the correlation coefficient between the cosine similarity and mean resident similarity on concept pairs
#from PMID: 16875881 
benchmark_similarity <- function(embedding,ref_cuis=NULL){
  results = data.frame(Similarity = 0,Relatedness = 0)
  #Take every CUI in the embedding if no reference specified 
  if(is.null(ref_cuis)){ref_cuis<-rownames(embedding)}
  #Load the file with concept pairs and mean resident scores 
  similar <- read.csv("./data/benchmarks/Similarity and Relatedness /UMNSRS_similarity.csv", header=TRUE)
  #Generate a data frame to store the cosine similarity and mean resident similarity for concept pairs
  df <- data.frame(mean=numeric(),cos=numeric())
  for(i in 1:length(similar$CUI1)){
    #Selecting only concept pairs where both CUIs are in the embedding
    if((similar$CUI1[i] %in% ref_cuis) & (similar$CUI2[i] %in% ref_cuis)){
      #Write to the data frame
      df[dim(df)[1]+1,]=c(similar$Mean[i],cos_similairty(embedding[similar$CUI1[i],],embedding[similar$CUI2[i],]))
    }
  }
  results$Similarity = cor(df$mean,df$cos,method="spearman")
  
  relatedness <- read.csv("./data/benchmarks/Similarity and Relatedness /UMNSRS_relatedness.csv", header=TRUE)
  #Generate a data frame to store the cosine similarity and mean resident similarity for concept pairs
  df <- data.frame(mean=numeric(),cos=numeric())
  for(i in 1:length(relatedness$CUI1)){
    #Selecting only concept pairs where both CUIs are in the embedding
    if((relatedness$CUI1[i] %in% ref_cuis) & (relatedness$CUI2[i] %in% ref_cuis)){
      #Write to the data frame
      df[dim(df)[1]+1,]=c(relatedness$Mean[i],cos_similairty(embedding[relatedness$CUI1[i],],embedding[relatedness$CUI2[i],]))
    }
  }
  results$Relatedness = cor(df$mean,df$cos,method="spearman")
  
  return(results)
}

