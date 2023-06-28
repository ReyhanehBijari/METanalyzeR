#' Rank CI Function
#'
#' Doing some analysis on the ranks and generating the desired output using mean phenotype method
#'
#' @param experiment_data input data
#' @param genotype_set vector of genotypes of interest
#' @param boots_matrix matrix of bootstrap resamples of locations (default = NULL).
#' By default it generates a bootstrap resample with 800 iteration.
#' @param p confidence interval (default = 0.8)
#' @param output_name Name of the file to be stored in /Outputs/ directory.
#' @import data.table
#' @import reshape2
#' @import dplyr
#' @import ggplot2
#' @import forcats
#' @import parallel
#' @export
#' @return A list of 6 items:
#' \itemize{
#'     \item rank_probabilities
#'     \item rank_info
#'     \item bs_ranks
#'     \item bs_PTs
#'     \item RI
#'     \item rankCI
#' }
#'

Rank_CI = function(experiment_data,
                         output_name = NULL,
                         genotype_set=NULL,
                         boots_matrix=NULL,
                         p = NULL){
  options(warn=-1)
  options(dplyr.summarise.inform = FALSE)
  suppressMessages({
    reshape2::melt(head(mtcars))
  })



  if (is.null(genotype_set)){
    genotype_set = unique(experiment_data$GENOTYPE)
  }

  if (is.null(p)){
    p = 0.8
  }

  if (is.null(output_name)){
    output_name = "output"
  }


  ##### Initial Preprocessing
  #experiment_data = experiment_data %>% filter(GENOTYPE %in% genotype_set)
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



    boots_list_obs = mclapply(1:nrow(boots_matrix), function(i){

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

    Rank_Interval_2 = function(CI_data ,
                               p
    ){


      #my_row = CI_data%>%dplyr::select(-c(GENOTYPE, Rank, M_Corrected_PT, MostProbable_rank))
      #my_row = my_row[2,]
      data_CI = CI_data[-which(names(CI_data)%in%c("GENOTYPE","M_Corrected_PT","Rank"))]
      RI= t(apply(data_CI,
                  1 ,
                  function(my_row){
                    #print("---------------------------")

                    non_zero = which(my_row !=0)
                    avail_prob = as.numeric(my_row[non_zero])
                    max_prob = max(which(avail_prob == max(avail_prob)))


                    left  = T
                    right = T

                    sums = unique(avail_prob[max_prob])
                    w = max_prob
                    id = max_prob
                    avail_prob[id] = -1


                    while(sums < p & ( left == T  |right == T)){

                      left  = ifelse(left  == T & id-1 >=1,T,F)
                      right = ifelse(right == T & id+1 <= length(avail_prob),T,F)

                      pool = which(avail_prob>0)

                      if(left == T){
                        L_id = pool[max(which(pool < id))]
                        L_prob = avail_prob[L_id]
                      }else{
                        L_prob = -1
                      }

                      if(right == T){
                        R_id = pool[min(which(pool > id))]
                        R_prob = avail_prob[R_id]
                      }else{
                        R_prob = -1
                      }

                      if(is.na(R_prob) == T | is.na(L_prob) == T) break
                      if(R_prob > L_prob){
                        sums = sums + R_prob
                        avail_prob[R_id] = -1
                        w = c(w , R_id)
                        id = R_id
                      }
                      if(R_prob < L_prob){
                        sums = sums + L_prob
                        avail_prob[L_id] = -1
                        w = c(w , L_id)
                        id = L_id
                      }
                      if(R_prob == L_prob){
                        if(R_prob > 0 ){
                          sums = sums + R_prob
                          avail_prob[R_id] = -1
                          w = c(w , R_id)
                          id = R_id
                        }
                      }

                    }#end while

                    ans = my_row
                    min_r = min(non_zero[w])
                    max_r = max(non_zero[w])
                    #print(paste("min:",min_r,"max:",max_r))

                    ans[-(min_r:max_r)]=NA
                    ans
                  }))
      RI = as.data.frame(RI)
      return(RI)


    }

    #################################
    #RI = Rank_Interval_1(CI_data,p)
    RI = Rank_Interval_2(CI_data, p)

    RI$GENOTYPE = CI_Info$GENOTYPE

    gg_df = melt(RI,
                 id.vars = "GENOTYPE" ,
                 variable.name = "RANK" ,
                 value.name = "PROB"

    )

    gg_df$RANK = as.numeric(as.character(gg_df$RANK))

    gg_df = gg_df %>% merge(CI_Info , by = "GENOTYPE")


    gg_df11 = merge(aggregate(PROB~GENOTYPE, gg_df, max ),gg_df)%>%unique()
    gg_df111 =gg_df11[-3]%>%unique()

    # sorting type is MostProbable_rank
    gg_df111 = gg_df111[order(gg_df111$MostProbable_rank, -gg_df111$PROB),]

    Var_Level = unique(gg_df111$GENOTYPE)

    ## MAIN SORTING
    gg_df$GENOTYPE = factor(gg_df$GENOTYPE , levels = rev(Var_Level))

    gg_df = gg_df[which(gg_df$MostProbable_rank <= length(unique(experiment_data$GENOTYPE))),]

    gg_df = gg_df[order(gg_df$MostProbable_rank),]


    # Plot
    gg = ggplot() +
      geom_tile(data = gg_df,
                aes(x=RANK,
                    y = GENOTYPE,
                    fill = -PROB)
      )+
      theme(panel.border = element_rect(colour = "black",size = 1,
                                        linetype = 1, fill = NA))+ #axis.text.x = element_text(angle = 90),
      scale_x_continuous(expand = c(0, 0))+

      geom_text(data = gg_df ,
                aes(x=RANK , y = GENOTYPE,label=PROB, fontface = 'bold'),
                angle = 0,
                hjust=0.5,colour="white",size=4
      )+

      scale_fill_gradient(guide=F,na.value = "gray93")+
      labs( y = "Genotypes")+ #title = paste0(p*100,"% Confidence Interval - ", exp),
      theme(
        plot.title = element_text(
          size = 10,
          face = "bold",
          hjust = 0.5
        ),axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))#+coord_equal()
    gg


rankCI = gg
rank_probabilities = CI_data[,!names(CI_data) %in% c("M_Corrected_PT", "Rank")]
rank_info = cbind(CI_Info,CI_data[,c("M_Corrected_PT","Rank")])
rank_info <- rank_info[ , c("GENOTYPE", "M_Corrected_PT", "Rank", "MostProbable_rank")]





  # Check if the directory exist:
  if (!file.exists(paste0('Outputs/', output_name, '.Rdata'))) {
    cat("A new folder is created as /Outputs to store the outputs.")
    dir.create('Outputs', showWarnings = FALSE)
  }






output = list(experiment_data = experiment_data,
                  rank_probabilities = rank_probabilities ,
                  rank_info = rank_info,
                  bs_ranks = boost_df_obs,
                  bs_PTs = boost_PT_info_obs,
                  RI = RI,
                  rankCI = gg)



save(output, file=paste0('Outputs/', output_name, '.Rdata'))

print(rankCI)


return(output)

}
