#' Rank Analyzer Function
#'
#' Doing some analysis on the ranks and generating the desired output
#'
#' @param experiment_data input data
#' @param genotype_set vector of genotypes of interest
#' @param boots_matrix matrix of bootstrap resamples of locations (default=NULL).
#' By default it generates a bootstrap resample with 1000 iteration.
#' @param method the method to do the analysis based on it
#' @param output_name Name of the file to be stored in /Outputs/ directory.
#' @import data.table
#' @import reshape2
#' @import dplyr
#' @export
#' @return A list of 6 items:
#' \itemize{
#'     \item CI_data
#'     \item CI_Info
#'     \item boost_df_obs
#'     \item boost_PT_info_obs
#'     \item pairwise_probabilities
#'     \item probabilistic_ranks
#' }
#'

rank_analyzer = function(experiment_data,
                         output_name,
                         genotype_set="",
                         boots_matrix=NULL,
                         method='both'){

  if (!(method %in% c('probabilistic', 'mean_phenotype', 'both'))){
    stop("Wrong Input: method should be one of the probabilistic, mean_phenotype, or both")
  }


  if (is.null(boots_matrix)){
    seed = 20
    iteration = 10

    set.seed(seed)
    N_env = length(unique(experiment_data$ENVIRONMENT))
    boots_matrix = lapply(1:iteration, function(i){
      as.character(sample(N_env,
                          N_env,
                          replace = T))
    }
    )

    index_out = do.call("rbind", boots_matrix)

    boots_matrix = index_out[1:iteration,]
  }






  ##### Initial Preprocessing
  include_v = unique(experiment_data$GENOTYPE)
  #if(checks == F){include_v = include_v[which(!include_v %in% is_check)]}

  print(paste0("included genotypes=", length(unique(include_v))))

  N_var = length(unique(experiment_data$GENOTYPE))
  print(paste0("N_var=",N_var))

  N_env = ncol(boots_matrix)

  INFO = RankPhenotype_extractor(experiment_data, unique(experiment_data$ENVIRONMENT))

  Rank_Table  = INFO[[1]]
  for(i in 2:ncol(Rank_Table)){
    Rank_Table[,i]=as.numeric(Rank_Table[,i])
  }

  Rank_Table=Rank_Table[(order(Rank_Table$GroundTruth)), ]

  PT_Table = INFO[[2]]



  if ((method == 'probabilistic') | (method == 'both')){


    boots_train_prob = list()

    boots_list_prob = lapply(1:nrow(boots_matrix), function(i){
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
        it_train_m = rbindlist(it_train_list)
        it_train_m

      })
      boots_train_prob = rbindlist(boots_train_la)
      boots_train_prob
    })
    boots_out_prob = data.table::rbindlist(boots_list_prob)

    pairwise_probabilities = boots_out_prob%>%
      group_by(GENOTYPE,Compared_to)%>%
      dplyr::summarise(Probability = round(sum(Prob,na.rm=T)/(N_env*nrow(boots_matrix)),3))

    pairwise_probabilities$Compared_to = as.character(pairwise_probabilities$Compared_to)

    probabilistic_ranks = Probabilistic_ranks(experiment_data, pairwise_probabilities)

  }







  if((method == 'mean_phenotype') | (method == 'both')){

    boots_list_obs = lapply(1:nrow(boots_matrix), function(i){

      tr_index_obs =  boots_matrix[i,]
      #tr_index_obs =  as.character(sample(N_env, N_env, replace = T))


      rank_table_train_obs = as.data.frame(Rank_Table$GENOTYPE)
      PT_table_train_obs = as.data.frame(Rank_Table$GENOTYPE)
      names(PT_table_train_obs)[1] = "GENOTYPE"

      for(k in 1:N_env){
        rank_table_train_obs = cbind(rank_table_train_obs,
                                     Rank_Table %>% dplyr::select(
                                       paste0("Rank_",tr_index_obs[k]))
        )

        PT_table_train_obs = cbind(PT_table_train_obs,
                                      PT_Table %>% dplyr::select(
                                        paste0("MeanPT_",tr_index_obs[k]))
        )
      }

      names(rank_table_train_obs) = c("GENOTYPE", tr_index_obs)

      if (!('EXPERIMENT' %in% names(experiment_data))){
        experiment_data['EXPERIMENT'] = 'EXPERIMENT1'
      }

      Experiments = unique(experiment_data$EXPERIMENT)
      Experiment = unique(experiment_data$EXPERIMENT)

      it_train_obs   = gen_comparison(experiment_data, Experiments, rank_table_train_obs)[[2]]


      PT_info_obs = data.frame("GENOTYPE" = PT_table_train_obs$GENOTYPE,
                                  "MeanPT" = round(rowSums(
                                    PT_table_train_obs[,-1])/N_env,3))

      #boots_train [[i]]= it_train_obs
      #tr_index1 [[i]]= data.frame(tr_index_obs)
      list(it_train_obs, PT_info_obs)
    })

    boost_df_obs = NULL
    boost_PT_info_obs = NULL

    for(i in 1:nrow(boots_matrix)){
      boost_df11_obs=boots_list_obs[[i]][[1]]
      boost_PT_info11_obs =boots_list_obs[[i]][[2]]
      boost_df_obs = rbind(boost_df_obs,boost_df11_obs)
      boost_PT_info_obs = rbind(boost_PT_info_obs,boost_PT_info11_obs)
    }


    colnames(boost_df_obs)[1] = "GENOTYPE"


    #currently excluding checks
    prob_df_obs = genotype_list_probs(experiment_data,boost_df_obs)
    CI_data = CI_calculator(experiment_data,prob_df_obs,perc)
    CI_Info = CI_data %>% dplyr::select(c(GENOTYPE,MostProbable_rank))
    order_obs = CI_Info$GENOTYPE

    CI_Info = CI_Info%>%merge(experiment_data[,c("GENOTYPE")], by= "GENOTYPE")%>%unique()
    CI_Info = CI_Info[match(order_obs,CI_Info$GENOTYPE), ]
    CI_data = CI_data %>% dplyr::select(-c(MostProbable_rank))

  }


  # Check if the directory exist:
  if (!file.exists(paste0('Outputs/', output_name, '.Rdata'))) {
    cat("A new folder is created as /Outputs to store the outputs.")
    dir.create('Outputs', showWarnings = FALSE)
  }






  if (method == 'probabilistic'){
    output = list(experiment_data = experiment_data,
                  pairwise_probs = pairwise_probabilities,
                  probabilistic_ranks = probabilistic_ranks)
  } else if(method == 'mean_phenotype'){
    output = list(experiment_data = experiment_data,
                  CI_data = CI_data ,
                  CI_Info = CI_Info,
                  bs_ranks = boost_df_obs,
                  bs_PTs = boost_PT_info_obs)
  } else {
    output = list(experiment_data = experiment_data,
                  CI_data = CI_data,
                  CI_Info = CI_Info,
                  bs_ranks = boost_df_obs,
                  bs_PTs = boost_PT_info_obs,
                  pairwise_probs = pairwise_probabilities,
                  probabilistic_ranks = probabilistic_ranks)
  }


  save(output, file=paste0('Outputs/', output_name, '.Rdata'))
  return(output)

}
