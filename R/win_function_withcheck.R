#' Win function with check
#'
#' do the same as gene_comparison but with checks
#'
#' @param Rank_Table_train1 Input data
#' @param var1 variable 1
#' @import dplyr
#' @import reshape2
#'
#' @export
#' @return pdata
#'
win_function_withcheck = function(Rank_Table_train1,
                                  var1
){

  var_to_analyze = setdiff(Rank_Table_train1$GENOTYPE,var1)
  wins = data.frame(matrix(0, nrow = 1,ncol = ncol(Rank_Table_train1)))
  names(wins) = colnames(Rank_Table_train1)
  wins$GENOTYPE = var_to_analyze

  #Computing Win Probability of GENOTYPE in each environment

  Rank_Table_train1[,-1] = apply(Rank_Table_train1[,-1],c(1,2) , as.character)
  Rank_Table_train1[,-1] = apply(Rank_Table_train1[,-1],c(1,2) , as.numeric)

  #win_data1 = data.frame("GENOTYPE" = Rank_Table_train1$GENOTYPE[which(Rank_Table_train1$GENOTYPE=="V143929")],
  #                      "BETTERTHAN" = Rank_Table_train1$GENOTYPE[which(Rank_Table_train1$GENOTYPE==var)],
  #                      "Prob" = 0)


  comp_indexx = which(Rank_Table_train1$GENOTYPE%in%var1)
  toanalyze_indexx = setdiff(1:2,comp_indexx)
  #print(comp_indexx)
  pdata = as.data.frame(sapply(2:ncol(wins), function(x) {

    wins[,x] =

      ifelse(Rank_Table_train1[comp_indexx,x] >=

               Rank_Table_train1[toanalyze_indexx,x],1,0)

  })
  )
  # for(i in 1:length(Ranks_table$GENOTYPE)){

  names(pdata) = var_to_analyze

  return(pdata)

}
