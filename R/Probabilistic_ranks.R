#' Rank Analyzer Function
#'
#' Doing some analysis on the ranks and generating the desired output
#'
#' @param subset_exp input data
#' @param pairwise_probabilities pairwise comparison probabilities
#' @import dplyr
#' @import reshape2
#' @export
#' @return A table of probabilistic ranks
#'
Probabilistic_ranks= function(subset_exp,
                              pairwise_probabilities

){

  realdata_d = subset_exp %>%
    group_by(GENOTYPE) %>%
    mutate(Mean_ObservedPT = round(mean(PT),3))%>%
    dplyr::select("GENOTYPE",
                  "Mean_ObservedPT")%>% unique()

  realdata_d$observed_Rank = rank(-round(realdata_d$Mean_ObservedPT,3),
                                  ties.method = "min")

  pairwise_probabilities = pairwise_probabilities %>% merge(realdata_d,by="GENOTYPE")

  pairwise_probabilities$ProbabilisticRank = pairwise_probabilities$observed_Rank


  iter_data = pairwise_probabilities

  if(length(unique(iter_data$ProbabilisticRank))!=length(unique(iter_data$GENOTYPE))){
    a = setdiff(1:length(unique(iter_data$GENOTYPE)), unique(iter_data$ProbabilisticRank))

    for(ii in a){
      print(ii)
      if(length(setdiff(1:length(unique(iter_data$GENOTYPE)),unique(iter_data$ProbabilisticRank)))==0){
        break
      }
      code = unique(iter_data$GENOTYPE[which(iter_data$ProbabilisticRank==(ii-1))])

      if(length(code)==0){
        u1 = 0
      }else if(length(code)==2){


        if(iter_data$Probability[which(iter_data$GENOTYPE==code[1]&iter_data$Compared_to==code[2])]>0.5){
          iter_data$ProbabilisticRank[which(iter_data$GENOTYPE==code[2])]=ii
        }else{
          iter_data$ProbabilisticRank[which(iter_data$GENOTYPE==code[1])]=ii
        }


      }else {
        num = length(code)
        r = 0
        for(j in 1:num){

          iter_data$ProbabilisticRank[which(iter_data$GENOTYPE==code[j])]=  ii-1+r
          print(ii-1+r)
          r = r+1
        }


        swaps=1
        loop=0
        while(swaps!=0){
          loop=loop+1
          #print("1")
          swaps = 0
          #print("2")

          for(k in (ii-1):(ii+num-1-2)){
            #print("3")

            compvar = unique(iter_data$GENOTYPE[which(iter_data$ProbabilisticRank==k)])
            #print("4")

            var = unique(iter_data$GENOTYPE[which(iter_data$ProbabilisticRank==k+1)])
            # print("5")

            if(iter_data$Probability[which(iter_data$ProbabilisticRank==k+1&iter_data$Compared_to==compvar)]>0.5){
              # print("6")

              iter_data$ProbabilisticRank[which(iter_data$ProbabilisticRank==k)]=k+1
              #print("7")

              iter_data$ProbabilisticRank[which(iter_data$ProbabilisticRank==k+1&iter_data$GENOTYPE==var)]=k
              #print("8")

              swaps=swaps+1
              print(swaps)

            }

          }

        }


      }
    }
  }


  iter_data = iter_data[order(iter_data$ProbabilisticRank),]

  swaps=1
  loop=0
  while(swaps!=0){
    loop=loop+1
    #print("1")
    swaps = 0
    #print("2")

    for(k in 1:(length(unique(iter_data$GENOTYPE))-1)){
      #print("3")

      compvar = unique(iter_data$GENOTYPE[which(iter_data$ProbabilisticRank==k)])
      #print("4")

      var = unique(iter_data$GENOTYPE[which(iter_data$ProbabilisticRank==k+1)])
      # print("5")

      if(iter_data$Probability[which(iter_data$ProbabilisticRank==k+1&iter_data$Compared_to==compvar)]>0.5){
        # print("6")

        iter_data$ProbabilisticRank[which(iter_data$ProbabilisticRank==k)]=k+1
        #print("7")

        iter_data$ProbabilisticRank[which(iter_data$ProbabilisticRank==k+1&iter_data$GENOTYPE==var)]=k
        #print("8")

        swaps=swaps+1
        print(swaps)

      }

    }

  }

  iter_data = iter_data %>% dplyr::select(c("GENOTYPE",
                                            "observed_Rank",
                                            "ProbabilisticRank")) %>% unique()


  iter_data = iter_data[order(iter_data$GENOTYPE),]
  iter_data$ProbabilisticRank = as.integer(iter_data$ProbabilisticRank)

  return(iter_data)

}
