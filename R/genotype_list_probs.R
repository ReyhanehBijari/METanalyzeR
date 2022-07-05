#' Distribution of rank probabilities
#'
#' Calculates how the rank probabilities are distributed
#'
#' @param subset_exp The input file
#' @param boots_train training bootstrap resamples
#' @import dplyr
#' @import reshape2
#' @export
#' @return probabilities for different varieties
#'
genotype_list_probs = function(subset_exp, boots_train)
  {
  #Var_nonchecks  = unique(subset_exp$GENOTYPE[which(subset_exp$CHECK=="FALSE")])
  Var_s  = unique(subset_exp$GENOTYPE)
  #print(length(Var_nonchecks))
  data = data.frame("GENOTYPE" = as.character(),
                    "Rank" = as.numeric() ,
                    "Prob" = as.numeric()
  )
  Nonchecks = FALSE

  myprobs=lapply(1: length(Var_s), function(u){
    var_set = boots_train %>% filter(GENOTYPE %in% Var_s[u])

    #prob_f should be changed to changed to prob_f_noncheck for just non-check varieties


    x  = as.numeric(as.character(prob_f(subset_exp,var_set,u,Nonchecks)[[1]]))
    px = round(prob_f(subset_exp,var_set,u,Nonchecks)[[2]],4)


    #skew_var_set = skewness(prob_f(subset_exp,var_set,u)[[2]])
    Ex_var_set = sum(x * px)
    Vx_var_set = sum(((x-Ex_var_set)^2)*px)

    probs = data.frame("GENOTYPE" = Var_s[u],
                       "Rank"    = x ,
                       "Prob"    = px
    )
    probs
  }
  )

  data = as.data.frame(rbindlist(myprobs))

  var_probs = dcast(data,
                    GENOTYPE~Rank,
                    fun.aggregate = sum
  )

  return(var_probs)

}
