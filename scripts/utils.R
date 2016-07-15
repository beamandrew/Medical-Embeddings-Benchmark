library(dplyr)
library(readr)
library(text2vec)
library(magrittr)
library(Rtsne)
library(ggplot2)
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

load_embeddings <- function(filename,convert_to_cui=TRUE,header=F,skip=1,sep=" ") {
  embeddings <- read.delim(filename,sep=sep,skip=skip,header=header)
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

load_semantic_type <- function(filename) {
  semantic <- read.delim(filename)
  return(semantic)
}

load_causitive <- function(filename){
  causitive <- read.delim(filename)
  return(causitive)
}

load_ndf_rt <- function(filename){
  ndf_rt <- read.delim(filename)
  return(ndf_rt)
}


benchmark_map <- function(embedding,k,ref_embeddings=NULL,take_intersection=TRUE){
  df <- data.frame(test = character(), embedding_name = character(), score = numeric(), stringsAsFactors = FALSE)
  ref_cuis <- rownames(embedding)
  if(!is.null(ref_embeddings)&take_intersection){
    for(j in 1:length(ref_embeddings)){
      ref_cuis <- intersect(ref_cuis,rownames(ref_embeddings[[j]]))
    }
  }
  #Benchmark the embedding 
  causitive <- benchmark_causitive(embedding,k,ref_cuis)
  semantic_type <- benchmark_semantic_type(embedding,k,ref_cuis)
  ndf_rt <- benchmark_ndf_rt(embedding,k,ref_cuis)
  
  name <- deparse(substitute(embedding))
  for(i in 1:dim(causitive)[1]){
    df[dim(df)[1]+1,] <- c(paste0('causitive_',causitive[i,1]),name,causitive[i,2])
  }
  for(i in 1:dim(semantic_type)[1]){
    df[dim(df)[1]+1,] <- c(paste0('semantic_type_',semantic_type[i,1]),name,semantic_type[i,2])
  }
  for(i in 1:dim(ndf_rt)[1]){
    df[dim(df)[1]+1,] <- c(paste0('ndf_rt_',ndf_rt[i,1]),name,ndf_rt[i,2])
  }
  
  #Benchmark the reference embeddings
  for(j in 1:length(ref_embeddings)){
    if(!take_intersection){
      ref_cuis <- intersect(rownames(embedding,rownames(ref_embeddings[[j]])))
    }
    causitive <- benchmark_causitive(ref_embeddings[[j]],k,ref_cuis)
    semantic_type <- benchmark_semantic_type(ref_embeddings[[j]],k,ref_cuis)
    ndf_rt <- benchmark_ndf_rt(ref_embeddings[[j]],k,ref_cuis)
    
    name <- paste0('reference_',j)
    for(i in 1:dim(causitive)[1]){
      df[dim(df)[1]+1,] <- c(paste0('causitive_',causitive[i,1]),name,causitive[i,2])
    }
    for(i in 1:dim(semantic_type)[1]){
      df[dim(df)[1]+1,] <- c(paste0('semantic_type_',semantic_type[i,1]),name,semantic_type[i,2])
    }
    for(i in 1:dim(ndf_rt)[1]){
      df[dim(df)[1]+1,] <- c(paste0('ndf_rt_',ndf_rt[i,1]),name,ndf_rt[i,2])
    }
  }
  
  df$score <- as.numeric(df$score)
return(df)
}

#Takes in one embedding you want to compare to some list of reference embeddings
benchmark_dcg<- function(embedding,k,ref_embeddings=NULL,take_intersection=TRUE,return_max=FALSE){
  df <- data.frame(test = character(), embedding_name = character(), score = numeric(), stringsAsFactors = FALSE)
  ref_cuis <- rownames(embedding)
  if(!is.null(ref_embeddings)&take_intersection){
    for(j in 1:length(ref_embeddings)){
      ref_cuis <- intersect(ref_cuis,rownames(ref_embeddings[[j]]))
    }
  }
  #Benchmark the embedding
  comorbidity <- benchmark_comorbidities(embedding, k, ref_cuis,return_max)
  name <- deparse(substitute(embedding))
  for(j in 1:dim(comorbidity)[1]){
      df[dim(df)[1]+1,] <- c(paste(comorbidity[j,1],comorbidity[j,2],sep='_'),name,comorbidity[j,3])
  }
  

  #Benchmark the reference embeddings
  for(j in 1:length(ref_embeddings)){
    if(!take_intersection){
      ref_cuis <- intersect(rownames(embedding,rownames(ref_embeddings[[j]])))
    }
    comorbidity <- benchmark_comorbidities(ref_embeddings[[j]], k, ref_cuis,return_max)
    name <- deparse(substitute(ref_embeddings[[j]]))
    for(i in 1:dim(comorbidity)[1]){
      df[dim(df)[1]+1,] <- c(paste(comorbidity[i,1],comorbidity[i,2],sep='_'),name,comorbidity[i,3])
    }
  }
  df$score <- as.numeric(df$score)
  return(df)
}


visualize_dcg <- function(df){
  rt <- ggplot(data=df, aes(x=test,y=score,color=embedding_name))+geom_point()
  rt <- rt+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  rt <- rt+labs(title='DCG Benchmark')+labs(x='Test')+labs(y='DCG Score')
  return(rt)
}

visualize_map <- function(df){
  rt <- ggplot(data=df, aes(x=df$test,y=df$score,color=df$embedding_name))+geom_point()
  rt <- rt+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  rt <- rt+labs(title='MAP Benchmark of')+labs(x='Test')+labs(y='Mean Average Precision')
  return(rt)
}

get_tsne <- function(embedding){
  return(Rtsne(as.matrix(embedding)))
}

visualize_embedding <- function(embedding,tsne=NULL,file='', type=''){
  type = toupper(type)
  if(is.null(tsne)){
    tsne <- get_tsne(embedding)
  }
  if(file==''){
    df <- tsne$Y
    colnames(df)<-c('X1','X2')
    rt <-ggplot(data=df,aes(x=X1,y=X2))+geom_point()+scale_alpha(.1)+theme_bw()
    return(rt)
  }
  name = strsplit(tail(strsplit(file,'/')[[1]],1),'.txt')[[1]]
  if(identical(type,'COMORBIDITY')){
    comorbidity <- load_comorbidity(file)
    df <- data.frame(tsne$Y)
    colnames(df)<-c('X1','X2')
    concepts <- intersect(comorbidity$CUI[which(comorbidity$Type=='Concept')],rownames(embedding))
    associations <- intersect(comorbidity$CUI[which(comorbidity$Type=='Association')],rownames(embedding))
    df$CUI <- rownames(embedding)
    df$Type <- 'Background'
    for(cui in concepts){
      df[match(cui,df$CUI),4]<-'Concept'}
    for(cui in associations){df[match(cui,df$CUI),4]<-'Association'}
    rt <- ggplot(data = df, aes(x=X1,y=X2,color=Type,alpha=Type))+geom_point()+scale_color_manual(values=c("#E69F00","#999999", "#56B4E9"))+scale_alpha_manual(values = c(1,.1,1))+theme_bw()
    rt <- rt + labs(title=paste(name,'Comorbidity t-SNE'))
    return(rt)
  }
  if(identical(type,'CAUSATIVE')){
    cause <- load_causitive(file)
    df <- data.frame(tsne$Y)
    colnames(df)<-c('X1','X2')
    df$CUI <- rownames(embedding)
    df$Type <- 'Background'
    for(cui in intersect(cause$CUI_Cause,df$CUI)){df[match(cui,df$CUI),4]<-'Cause'}
    for(cui in intersect(cause$CUI_Result,df$CUI)){df[match(cui,df$CUI),4]<-'Result'}
    rt <- ggplot(data = df, aes(x=X1,y=X2,color=Type,alpha=Type))+geom_point()+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+scale_alpha_manual(values = c(.1,1,1))+theme_bw()
    rt <- rt+labs(title=paste(name,'Causitive t-SNE'))
    return(rt)
  }
  if(identical(type,'SEMANTIC_TYPE')){
    semantic_type <- load_semnatic_type(file)
    df <- data.frame(tsne$Y)
    colnames(df)<-c('X1','X2')
    cuis <- intersect(semantic_type$CUI,rownames(embedding))
    df$CUI <- rownames(embedding)
    df$Type <- 'Background'
    for(cui in cuis){df[match(cui,df$CUI),4]<-name}
    rt <- ggplot(data = df, aes(x=X1,y=X2,color=Type,alpha=Type))+geom_point()+scale_color_manual(values=c("#999999", "#E69F00"))+scale_alpha_manual(values = c(.1,1))+theme_bw()
    rt <- rt+labs(title=paste(name,'Semantic Type t-SNE'))
    return(rt)
  }
  if(identical(type,'NDF_RT')){
    ndfrt <- load_semnatic_type(file)
    treatment <- intersect(ndfrt$Treatment,rownames(embedding))
    condition <- ''
    for (i in 1:length(ndfrt$Condition)){
      condition <- paste0(condition, ndfrt$Condition[i])
    }
    condition <- intersect(strsplit(condition,','),rownames(embedding))
    df <- data.frame(tsne$Y)
    colnames(df)<-c('X1','X2')
    df$CUI <- rownames(embedding)
    df$Type <- 'Background'
    for(cui in treatment){
      df[match(cui,df$CUI),4]<-'Treatment'}
    for(cui in condition){
      df[match(cui,df$CUI),4]<-'Condition'}
    rt <- ggplot(data = df, aes(x=X1,y=X2,color=Type,alpha=Type))+geom_point()+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+scale_alpha_manual(values = c(.1,1,1))+theme_bw()
    rt <- rt+labs(title=paste(name,'NDF_RT t-SNE'))
    return(rt)
  }
  #No valid type
  return(NULL)
}





