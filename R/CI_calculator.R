#' Confidence Interval calculator
#'
#' Calculate the CI for the probabilities
#'
#' @param exp_data Input data
#' @param experiment_data experiment data
#' @param p A float between 0 and 1 indicating the percentage used for constructing the CI
#' @import dplyr
#' @import reshape2
#' @export
#' @return The CI for the probabilities


CI_calculator = function(exp_data, experiment_data, p)
  {
  exp_data =  exp_data %>%
    group_by(GENOTYPE) %>%
    summarize(M_Corrected_PT = mean(M_Corrected_PT))


  experiment_data = merge(experiment_data, exp_data, by = "GENOTYPE")

  experiment_data$Rank = rank(-round(experiment_data$M_Corrected_PT,2),
                              ties.method = "min"
  )  #CHANGE

  experiment_data = experiment_data %>%
    mutate_if(is.numeric, round, digits = 2) #rounding off

  ee = ncol(experiment_data)
  #setting it to be the one with highest prob
  iddd = which(names(experiment_data)=="M_Corrected_PT")
  max_id = max.col(experiment_data[2:(iddd-1)],ties.method="first")
  experiment_data$MostProbable_rank  = colnames(experiment_data[2:(iddd-1)])[max_id]

  experiment_data$MostProbable_rank = gsub("X","",experiment_data$MostProbable_rank )

  experiment_data$MostProbable_rank = as.numeric(experiment_data$MostProbable_rank)
  experiment_data = experiment_data[order(experiment_data$MostProbable_rank),]

  return(experiment_data)


  ### Creating rank intervals #############
  #  rank_intervals = t(apply(experiment_data2 , 1 , function(my_row){

  #    ord = order(-my_row )             # order of the elements in the row
  #    ord_row = as.numeric(my_row[ord])
  #    sums = cumsum(ord_row)            # cumulative sumation by order
  #    w = min(which(sums >=p))       # where the cummulative sum passes p
  #    ans = my_row                      # a copy of row
  #    ans[-ord[1:w]]=0                  # replacing elements aren't used to zero
  #    ans
  #  }))

  #  return(rank_intervals)

}
