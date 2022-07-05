#' Pairwise Comparison
#'
#' Creates pairwise probability comparison
#'
#' @param subset_exp Input data
#' @param Rank_Table table of ranks
#' @param index_file matrix of the bootstrap resample of location
#' @param genotype_set vector of genotypes of interest
#' @export
#' @import dplyr
#' @import reshape2
#'
#'
#' @return pairwise comparison
#'
pairwise_probs = function(subset_exp,
                           Rank_Table,
                           index_file,
                           genotype_set
){

    N_env = ncol(index_file)
    #tr_index1 = list()
    boots_train = list()
    boots_train_result = NULL
    boots_train_result1 = NULL
    boots_train_result_total = NULL

    boots_list = lapply(1:nrow(index_file), function(i){
      tr_index = index_file[i,]
      Rank_Table_train = as.data.frame(Rank_Table$VARIETY)
      for(k in 1:N_env){
        Rank_Table_train = cbind(Rank_Table_train,
                                 Rank_Table%>%select(paste0("Rank_",tr_index[k]))
        )
      }
      names(Rank_Table_train) = c("VARIETY", tr_index)
      #Experiments = Locs

      boots_train_la = lapply(1: (length(genotype_set)), function(kk) {
        var_z = genotype_set[kk]
        it_train_list = lapply(setdiff(1: (length(genotype_set)),kk), function(z){  #setdiff(1: (length(genotype_set)-1)
          to_keep = which(Rank_Table_train$VARIETY%in%genotype_set[kk]|Rank_Table_train$VARIETY%in%genotype_set[z])
          Rank_Table_train1 = Rank_Table_train[to_keep,]
          it_train   = win_function_withcheck(Rank_Table_train1,genotype_set[z])
          it_train$VARIETY = names(it_train)[1]
          names(it_train)[1]="Prob"
          it_train$Compared_Var = genotype_set[z]
          it_train

          sum_prob = sum(it_train$Prob)
          out = it_train[1,]
          out$Prob = sum_prob
          out

        })
        it_train_m = rbindlist(it_train_list)
        it_train_m

      })
      boots_train = rbindlist(boots_train_la)
      boots_train
    })
    boots_out = rbindlist(boots_list)

    boots_train_result_total = boots_out%>%
      group_by(VARIETY,Compared_Var)%>%
      dplyr::summarise(Prob = sum(Prob,na.rm=T)/(N_env*nrow(index_file)))


    #colnames(boots_train)[1] = "VARIETY"
    return(boots_train_result_total)


  }
