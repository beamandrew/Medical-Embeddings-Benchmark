library(readr)
library(text2vec)
library(dplyr)
source("./Scripts/utils.R")
options(stringsAsFactors = FALSE)

benchmark_comorbidities <- function(embeddings,k) {
  df <- data.frame(class = character(), concept = character(), dcg = numeric(), stringsAsFactors = FALSE)
  for(file in list.files('./data/benchmarks/comorbidities/')){
    comorbidity <- load_comorbidity(paste0('./data/benchmarks/comorbidities/',file))
    concepts = which(comorbidity$Type=='Concept')
    for(i in 1:length(concepts)){
      print(comorbidity$CUI[concepts[i]])
      sim_scores <- get_dist(embeddings,comorbidity$CUI[concepts[i]])
      df[i,]<-c(strsplit(file,'.txt'), comorbidity$String[concepts[i]], dcg(sim_scores[1:k],comorbidity$CUI))
    }
  }
  return(df)
}