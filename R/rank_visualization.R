#' Rank Visualizer
#'
#' Plots the heatmap of the CI's for probabilities
#'
#' @param input_name Name of the saved file from `rank_analyzer` function which is stored
#' @param p A float in [0,1] used for constructing the confidence interval
#' @param n_top An integer value showing the number of top genotypes to be displayed (default=2)
#' in /Outputs/ directory.
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import forcats
#' @export
#' @return a ggplot object of the heatmap for the CI of the ranks.
#'
rank_visualization = function(input_name,
                              p,
                              n_top=2,
                              method='mean_phenotype'
){

  # Load CI data and CI Info


  if (method == 'mean_phenotype'){
    if (!file.exists(paste0('Outputs/', input_name, '.Rdata'))) {
      cat("The file doesn't exist. First run the `rank_analyzer` function to generate the input for this function.")
    } else {
      load(paste0('Outputs/', input_name, '.Rdata'))
      CI_data = output['CI_data']
      CI_Info = output['CI_Info']
      bs_ranks = output['bs_ranks']
      bs_PTs = output['bs_PTs']
      pairwise_probs = output['pairwise_probs']
      probabilistic_ranks = output['probabilistic_ranks']

    }

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
    #################################

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

    gg_df = gg_df[which(gg_df$MostProbable_rank <= n_top),]

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



    return(gg)

  } else if (method == 'probabilistic'){

      N_env = length(unique(experiment_data$ENVIRONMENT))
      original_environments = unique(experiment_data$ENVIRONMENT)
      experiment_data_backup = experiment_data

      seed = 20
      set.seed(seed)

      iteration_for_prob = 10

      index_list = lapply(1:25000, function(i){
        as.character(sample(N_env,
                            N_env,
                            replace = T))
      })

      index_out = do.call("rbind", index_list)

      index_file = index_out[1:iteration_for_prob,]



      myfile = NULL

      boots_probabilistic_list = lapply(1:nrow(index_file), function(i){

        sample_location_size=1
        loc_index = as.numeric(index_file[indexx,])


        include_env = original_environments[loc_index]
        include_env

        ex_dlist = lapply(loc_index, function(j){
          rand = sample(1:8000,1)
          ex_d = experiment_data%>%dplyr::filter(ENVIRONMENT==original_environments[j])
          ex_d$ENVIRONMENT = paste(ex_d$ENVIRONMENT,rand)
          ex_d
        })

        subset_exp = rbindlist(ex_dlist)

        unique(subset_exp$ENVIRONMENT)
        length(unique(subset_exp$ENVIRONMENT))

        subset_exp = subset_exp %>%
          group_by(GENOTYPE) %>% unique() %>%
          dplyr::mutate(N_ENV = length(unique(ENVIRONMENT)))


        length(main_v)

        length(unique(subset_exp$ENVIRONMENT))

        subset_exp = subset_exp%>%dplyr::filter(GENOTYPE%in%main_v)


        env_list = lapply(main_v, function(v){
          unique(subset_exp$ENVIRONMENT[which(subset_exp$GENOTYPE == v)])
        })
        length(Reduce(intersect, env_list))


        envv = Reduce(intersect, env_list)

        length(unique(subset_exp$GENOTYPE))

        subset_exp = subset_exp%>%dplyr::filter(ENVIRONMENT%in%envv)

        length(unique(subset_exp$GENOTYPE))

        Target_exp = unique(subset_exp$EXPERIMENT)



        #Number of locations to be observed
        N_env = length(unique(subset_exp$ENVIRONMENT))* sample_location_size
        Locs = unique(subset_exp$ENVIRONMENT)

        set.seed(seed)

        loc_index = sample(length(unique(subset_exp$ENVIRONMENT)),N_env, replace = F)
        loc_index

        include_env = unique(subset_exp$ENVIRONMENT)[loc_index]
        include_env

        ex_dlist = lapply(1:length(Locs), function(j){
          ex_d = subset_exp%>%dplyr::filter(ENVIRONMENT==Locs[j])
          ex_d
        })

        exp_data = rbindlist(ex_dlist)
        exp_data$GENOTYPE = as.factor(exp_data$GENOTYPE)



        genotype_set = unique(exp_data$GENOTYPE)#[1:10]

        exp_data_var = exp_data %>% dplyr::filter(GENOTYPE%in%genotype_set)
        exp_data = exp_data %>%
          group_by(GENOTYPE,ENVIRONMENT) %>%
          dplyr::mutate(REPNO = 1:length(ACTUAL_PT))

        Exp = unique(exp_data$EXPERIMENT)
        Target_exp = Exp


        p = 0.8
        CHECKS = "All"#"All" "Nonchecks"
        Nonchecks = FALSE  #If true, you are just considering nonchecks
        sorting = "MostProbable_rank" #Option: "blup", "MostProbable_rank"

        #for probability estimations
        included_env = "All"
        Output_Type  = "PT" #c("PT", "BLUP_OUTPUT")
        sample_location_size = 1



        ### data is the data before reducing the number of locations and/or
        #re sampling of locations
        ### exp_data is the data we are going to do bootstrapping on
        #***********************************************************************#

        ##### IMPORTANT: Individual Ranking Function is used

        realdata_d = exp_data %>%
          group_by(GENOTYPE) %>%
          mutate(Mean_ObservedPT = round(mean(ACTUAL_PT),3))
        realdata_d = realdata_d %>%
          group_by(GENOTYPE) %>%
          mutate(SD_ObservedYdiff = round(sd(ACTUAL_PT),3)) %>%
          dplyr::select(GENOTYPE,
                        Mean_ObservedPT,
                        SD_ObservedYdiff,
                        CHECK) %>%
          unique()

        realdata_d$Actual_Rank = rank(-round(realdata_d$Mean_ObservedPT,3),
                                      ties.method = "min")



        genotype_set = unique(exp_data$GENOTYPE)

        p = 0.8
        p = 0.8

        pairwise_probabilities = pairwise_probs(exp_data,
                                                checks = F,
                                                genotype_set,
                                                index_file)
        probabilistic_ranks = probabilistic_rank_analyzer(subset_exp, pairwise_probabilities)

        probabilistic_ranks

      })


      boots_probabilistic = rbindlist(boots_probabilistic_list)


      analys = boots_probabilistic %>% group_by(ProbabilisticRank,GENOTYPE) %>% mutate(sum_r = n())

      analys = analys %>% dplyr::select(GENOTYPE,ProbabilisticRank,sum_r)%>%unique()


      analys$PROB = round(analys$sum_r/nrow(mydata)*(length(unique(analys$GENOTYPE))),2)


      names(analys)[2] = "RANK"




      CI_MATRIX = as.data.frame(matrix(rep(0,length(unique(mydata$GENOTYPE))*length(unique(mydata$GENOTYPE))), nrow = length(unique(mydata$GENOTYPE)),  ncol = length(unique(mydata$GENOTYPE))))
      CI_MATRIX = cbind("GENOTYPE" = unique(analys$GENOTYPE),CI_MATRIX)



      for(k in 1:length(unique(analys$GENOTYPE))){
        for (i in unique(analys$RANK[which(analys$GENOTYPE==unique(CI_MATRIX$GENOTYPE)[k])])) {


          CI_MATRIX[k,i+1]= analys$PROB[which(analys$GENOTYPE==unique(CI_MATRIX$GENOTYPE)[k]&analys$RANK==i)]
        }
      }


      CI_data = CI_MATRIX
      p=0.8
      p=0.8
      n_top=length(unique(CI_data$GENOTYPE))
      sorting = "MostProbable_rank"

      Rank_Interval_2 = function(CI_data ,
                                 p
      ){

        #my_row = CI_data%>%dplyr::select(-c(GENOTYPE, Rank, M_Corrected_PT, MostProbable_rank))
        #my_row = my_row[2,]

        RI= t(apply(CI_data %>% dplyr::select(-1),
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
      #################################

      RI$GENOTYPE = CI_data$GENOTYPE
      names(RI)[-ncol(RI)] = 1:(ncol(RI)-1)
      gg_df = reshape2::melt(RI,
                             id.vars = "GENOTYPE" ,
                             variable.name = "RANK" ,
                             value.name = "PROB"

      )

      gg_df$RANK = as.numeric(as.character(gg_df$RANK))

      gg_df11 = merge(aggregate(PROB~GENOTYPE, gg_df, max ),gg_df)%>%unique()
      names(gg_df11)[3] = "MostProbable_rank"
      gg_df = gg_df %>% full_join(gg_df11[-2],by=c("GENOTYPE"))

      if(sorting == "MostProbable_rank") {

        gg_df11 = gg_df11[order(gg_df11$MostProbable_rank, -gg_df11$PROB),]

        Var_Level = unique(gg_df11$GENOTYPE)

        ## MAIN SORTING
        gg_df$GENOTYPE = factor(gg_df$GENOTYPE , levels = rev(Var_Level))

        gg_df = gg_df[which(gg_df$MostProbable_rank <= n_top),]

        gg_df = gg_df[order(gg_df$MostProbable_rank),]


      }


      gg = ggplot() +
        geom_tile(data = gg_df,
                  aes(x=RANK,
                      y = GENOTYPE,
                      fill = -PROB)
        )+
        theme(panel.border = element_rect(colour = "black",size = 1,
                                          linetype = 1, fill = NA))+ #axis.text.x = element_text(angle = 90),
        scale_x_continuous(expand = c(0, 0),breaks = 1:10)+

        geom_text(data = gg_df ,
                  aes(x=RANK , y = GENOTYPE,label=PROB, fontface = 'bold'),
                  angle = 0,
                  hjust=0.5,colour="white",size=4
        )+
        scale_fill_gradient(guide=F,na.value = "gray93")+
        labs( y = "Genotype")+ #title = paste0(p*100,"% Confidence Interval - ", exp),
        theme(
          plot.title = element_text(
            size = 10,
            face = "bold",
            hjust = 0.5
          ),axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
      gg
  }


}
