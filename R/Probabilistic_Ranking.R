#' Probabilistic Ranking Function
#'
#' Doing some analysis on the ranks and generating the desired output using the probabilistic method
#'
#' @param experiment_data input data
#' @param genotype_set vector of genotypes of interest
#' @param boots_matrix matrix of bootstrap resamples of locations (default = NULL).
#' By default it generates a bootstrap resample with 1000 iteration.
#' @param output_name Name of the file to be stored in /Outputs/ directory.
#' @import data.table
#' @import reshape2
#' @import dplyr
#' @export
#' @return A list of 2 items:
#' \itemize{
#'     \item pairwise_probabilities
#'     \item probabilistic_ranks
#' }
#'

Probabilistic_Ranking = function(experiment_data,
                         output_name,
                         genotype_set=NULL,
                         boots_matrix=NULL){
  options(warn=-1)
  options(dplyr.summarise.inform = FALSE)
  suppressMessages({
    reshape2::melt(head(mtcars))
  })



  if (is.null(genotype_set)){
    genotype_set = unique(experiment_data$GENOTYPE)
  }



  ##### Initial Preprocessing
  experiment_data = experiment_data %>% filter(GENOTYPE %in% genotype_set)
  include_v = unique(experiment_data$GENOTYPE)
  #if(checks == F){include_v = include_v[which(!include_v %in% is_check)]}

  print(paste0("included genotypes=", length(unique(include_v))))

  N_var = length(unique(experiment_data$GENOTYPE))
  print(paste0("N_var=",N_var))


  if (is.null(boots_matrix)){
    iteration = 800
    N_env = length(unique(experiment_data$ENVIRONMENT))
    boots_matrix = lapply(1:iteration, function(i){
      as.character(sample(N_env,
                          N_env,
                          replace = T))})

    index_out = do.call("rbind", boots_matrix)

    boots_matrix = index_out[1:iteration,]
  }


  # if (is.null(genotype_set)){
  #   genotype_set = unique(experiment_data$GENOTYPE)
  # }
  #
  #
  #
  # ##### Initial Preprocessing
  # experiment_data = experiment_data %>% filter(GENOTYPE %in% genotype_set)
  # include_v = unique(experiment_data$GENOTYPE)
  # #if(checks == F){include_v = include_v[which(!include_v %in% is_check)]}
  #
  # print(paste0("included genotypes=", length(unique(include_v))))
  #
  # N_var = length(unique(experiment_data$GENOTYPE))
  # print(paste0("N_var=",N_var))

  N_env = ncol(boots_matrix)

  INFO = RankPhenotype_extractor(experiment_data, unique(experiment_data$ENVIRONMENT))

  Rank_Table  = INFO[[1]]
  for(i in 2:ncol(Rank_Table)){
    Rank_Table[,i]=as.numeric(Rank_Table[,i])
  }

  Rank_Table=Rank_Table[(order(Rank_Table$GroundTruth)), ]

  PT_Table = INFO[[2]]


    boots_train_prob = list()

    boots_list_prob = mclapply(1:nrow(boots_matrix), function(i){
      tr_index_prob = boots_matrix[i,]
      Rank_Table_train_prob = as.data.frame(Rank_Table$GENOTYPE)
      for(k in 1:N_env){
        Rank_Table_train_prob = cbind(Rank_Table_train_prob,
                                      Rank_Table%>%dplyr::select(paste0("Rank_",tr_index_prob[k]))
        )
      }

      names(Rank_Table_train_prob) = c("GENOTYPE", tr_index_prob)
      #Experiments = Locs

      boots_train_la = lapply(1: (length(genotype_set)), function(kk) {
        var_z = genotype_set[kk]
        it_train_list = lapply(setdiff(1: (length(genotype_set)),kk), function(z){  #setdiff(1: (length(genotype_set)-1)
          to_keep = which(Rank_Table_train_prob$GENOTYPE%in%genotype_set[kk]|Rank_Table_train_prob$GENOTYPE%in%genotype_set[z])
          Rank_Table_train1 = Rank_Table_train_prob[to_keep,]
          it_train   = win_function_withcheck(Rank_Table_train1, genotype_set[z])
          it_train$GENOTYPE = names(it_train)[1]
          names(it_train)[1]="Prob"
          it_train$Compared_to = genotype_set[z]
          it_train

          sum_prob = sum(it_train$Prob)
          out = it_train[1,]
          out$Prob = sum_prob
          out

        })
        it_train_m = data.table::rbindlist(it_train_list)
        it_train_m

      })
      boots_train_prob = data.table::rbindlist(boots_train_la)
      boots_train_prob
    })
    boots_out_prob = data.table::rbindlist(boots_list_prob)

    pairwise_probabilities = boots_out_prob%>%
      group_by(GENOTYPE,Compared_to)%>%
      dplyr::summarise(Probability = round(sum(Prob,na.rm=T)/(N_env*nrow(boots_matrix)),3))

    pairwise_probabilities$Compared_to = as.character(pairwise_probabilities$Compared_to)

    probabilistic_ranks = Probabilistic_ranks(experiment_data, pairwise_probabilities)



  # Check if the directory exist:
  if (!file.exists(paste0('Outputs/', output_name, '.Rdata'))) {
    cat("A new folder is created as /Outputs to store the outputs.")
    dir.create('Outputs', showWarnings = FALSE)
  }



output = list(experiment_data = experiment_data,
                  pairwise_probs = pairwise_probabilities,
                  probabilistic_ranks = probabilistic_ranks)


save(output, file=paste0('Outputs/', output_name, '.Rdata'))


return(output)

}
