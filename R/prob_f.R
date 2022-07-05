#' probaility f
#'
#' helper function
#'
#' @param subset_exp inut data
#' @param new_analysis new analysis
#' @param u genotype to compare
#' @param Nonchecks Nonchecks
#' @import reshape2
#' @import dplyr
#' @export
#' @return rank & probs
#'

prob_f  = function(subset_exp,
                   new_analysis,
                   u,
                   Nonchecks = FALSE
){

  if(Nonchecks == TRUE){
    Var_s   = unique(subset_exp$GENOTYPE[which(subset_exp$CHECK=="FALSE")])
  }else
  {
    Var_s   = unique(subset_exp$GENOTYPE)
  }
  vars_ranks = new_analysis %>% dplyr::filter(
    GENOTYPE %in% Var_s[u])

  vars_ranks$group = as.character(c(1:nrow(vars_ranks)))

  vars_ranks_melted = reshape2::melt(
    vars_ranks[-1],
    id.GENOTYPE_list = "group",
    variable.name   = "Rank" ,
    value.name      = "Count"
  )

  vars_rankprobs = vars_ranks_melted %>%
    group_by(Rank) %>%
    dplyr::summarize(prob = sum(Count)/
                       (nrow(new_analysis)*
                          length(
                            unique(
                              subset_exp$ENVIRONMENT
                            )
                          )
                       )
    )

  vars_rankprobs$Rank = gsub("X","",vars_rankprobs$Rank)
  return(list(x  = vars_rankprobs$Rank,
              px = vars_rankprobs$prob
  )
  )
}
