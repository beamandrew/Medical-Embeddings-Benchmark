# takes in a list of strings that are the names of your embedding matrices in your local environment 

generate_BAT <- function(embedding_list,k,take_intersection=T){
  bat <- data.frame(Embedding=character(),Test=character(),Score=numeric(), stringsAsFactors = F)
  ref_cuis <- NULL
  #find the CUI intersection of all the embeddings
  if(take_intersection){
    for(i in 1:length(embedding_list)){
      if(i==1){
        ref_cuis <- rownames(get(embedding_list[[i]]))
      }
      ref_cuis <- intersect(rownames(get(embedding_list[[i]])),ref_cuis)
    }
  }
  print(paste('Total CUIs in intersection: ',length(ref_cuis)))
  for(i in 1:length(embedding_list)){
    embedding <- get(embedding_list[[i]])
    print(paste('Working on...',embedding_list[[i]]))
    #Benchmark comorbidities
    print('Benchmarking comorbidities...')
    df <- benchmark_comorbidities(embedding, k=k, ref_cuis, metric='AP', return_max = T)
    df$embedding <- embedding_list[[i]]
    bat <- combine_frames(bat,data.frame(cbind(df$embedding,df$class,df$score)))
    
    #Benchmark semantic type 
    print('Benchmarking semantic types...')
    df <- benchmark_semantic_type(embedding, k, ref_cuis)
    df$embedding <- embedding_list[[i]]
    bat <- combine_frames(bat, data.frame(cbind(df$embedding,df$semantic_type,df$map)))
    
    #Benchmark causitive 
    print('Benchmarking causitive...')
    df <- benchmark_causitive(embedding, k, ref_cuis)
    df$embedding <- embedding_list[[i]]
    bat <- combine_frames(bat, data.frame(cbind(df$embedding, df$cause, df$map)))
  
    #Benchmark NDF RT 
    print('Benchmarking NDF RT...')
    df <- benchmark_ndf_rt(embedding, k, ref_cuis)
    df$embedding <- embedding_list[[i]]
    bat <- combine_frames(bat, data.frame(cbind(df$embedding, df$ndf_rt, df$map)))
    
    #Benchmark similarity 
    print('Benchmarking similarity...')
    df <- benchmark_similarity(embedding, './data/benchmarks/Similarity and Relatedness /UMNSRS_similarity.csv', ref_cuis)
    bat[dim(bat)[1]+1,] <- c(embedding_list[[i]],'Similarity', df)
    
    #Benchmark relatedness
    print('Benchmarking relatedness...')
    df <- benchmark_similarity(embedding, './data/benchmarks/Similarity and Relatedness /UMNSRS_relatedness.csv', ref_cuis)
    bat[dim(bat)[1]+1,] <- c(embedding_list[[i]],'Relatedness', df)
    
    #Save after each embedding 
    save(bat,file='bat.Rda')
  }
  return(bat)
}

combine_frames <- function(master_frame, add_frame){
  if(is.null(master_frame)){
    return(add_frame)
  }
  if(is.null(add_frame)){
    return(master_frame)
  }
  for(i in 1:dim(add_frame)[1]){
    master_frame[dim(master_frame)[1]+1,] <- add_frame[i,]
  }
  return(master_frame)
}
