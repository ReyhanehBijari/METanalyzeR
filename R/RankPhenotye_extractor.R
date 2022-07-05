#' Rank & Phenotype Extractor
#'
#' Extracts the rank and yield tables
#'
#' @param subset_exp input data
#' @param include_env include environment or not
#' @import dplyr
#' @import reshape2
#' @export
#' @return Rank table & Yield table
#'
RankPhenotype_extractor = function (subset_exp,
                               include_env
){

  ### Calculating the ground truth

  gt =  subset_exp %>%
    group_by(GENOTYPE) %>%
    summarise(M_Corrected_PT = mean(M_Corrected_PT, na.rm = T))

  #Rank0 gives rank 1 to the min and largest for the max Ydiff
  gt$Rank = rank(-gt$M_Corrected_PT,
                 ties.method = "min") %>% as.factor()
  ground_truth_Rank = gt[,c(1,3)]
  names(ground_truth_Rank)[2] = "GroundTruth"

  ### Calculating Individual Ranking
  Rank_Table     = ground_truth_Rank
  #reordering the locations just in case!
  Environments = sample(include_env)
  PT_inf     = ground_truth_Rank[,1]

  for(n in 1:length(Environments)){
    env_data =  subset_exp %>%
      dplyr::filter(ENVIRONMENT == Environments[n]) %>%
      group_by(GENOTYPE) %>%
      dplyr::summarise(M_Corrected_PT = mean(Corrected_PT, na.rm = T))#!!!It is Corrected_PT that is changed to PT

    env_data$Rank = rank(-env_data$M_Corrected_PT,
                         ties.method = "min") %>% as.factor()

    names(env_data)[3] = paste0("Rank_",n)

    Rank_Table = merge(Rank_Table, env_data[,c(1,3)], all = TRUE)

    #Rank_Table  = Rank_Table %>% merge(env_data[,c(1,3)] , by = c("GENOTYPE"))
    names(env_data)[2] = paste0("MeanPT_",n)

    PT_inf = merge(PT_inf, env_data[,c(1,2)], all = TRUE)


  }
  #fct_explicit_na(Rank_Table, 0)
  Rank_Table = apply(Rank_Table,2, function(x) fct_explicit_na(x, na_level = "Missing"))
  Rank_Table = data.frame(Rank_Table)


  #fct_explicit_na(Rank_Table$Rank_2, na_level = "0")
  #Rank_Table[is.na(Rank_Table[,-c(1,2)])] <- 0
  return(list(Rank_Table,PT_inf))
}
