#' Probability Extractor
#'
#' Finds the probability of ranks for each genotype in all bootstrap resamples.
#'
#' @param Dataset Input data
#' @param Experiment The experiment
#' @param Ranks_table table of ranks
#' @import dplyr
#' @import reshape2
#' @export
#' @return a list of 2 items:
#' \itemize{
#'    \item df1
#'    \item Rey
#' }
#'
gen_comparison = function(Dataset,
                          Experiment,
                          Ranks_table
  ){

  wins = Ranks_table
  wins[,-1] = 0

  #Computing Win Probability of GENOTYPE in each environment

  Ranks_table[,-1] = apply(Ranks_table[,-1],c(1,2) , as.character)
  Ranks_table[,-1] = apply(Ranks_table[,-1],c(1,2) , as.numeric)

  for(j in 2:ncol(Ranks_table)){
    for(i in 1:length(Ranks_table$GENOTYPE)){

      wins[i,j] = length(which(Ranks_table[,j] > Ranks_table[i,j]))
    }
  }

  pdata = as.data.frame(t(wins[-1]))
  pdata = apply(pdata, 2, cumsum)


  for(i in 1:(nrow(pdata))){
    pdata[i,] = pdata[i,]/((nrow(Ranks_table)-1)*i)
  }

  pdata = as.data.frame(pdata)
  colnames(pdata) = as.character(Ranks_table$GENOTYPE)


  V_list = data.frame("GENOTYPE" = Ranks_table$GENOTYPE)

  dd =  Dataset %>%
    ungroup()%>%
    filter(GENOTYPE %in% V_list$GENOTYPE) %>%
    select(EXPERIMENT,GENOTYPE) %>%
    unique()

  dd =  V_list %>%
    merge(dd, by = "GENOTYPE")  %>%
    filter(EXPERIMENT==Experiment)


  #df = subset(pdata,dplyr::select = V_list$GENOTYPE)
  df = pdata

  df$NLoc = c(1:nrow(df))
  df1 = melt(df,
                       id.vars = "NLoc" ,
                       variable.name = "GENOTYPE" ,
                       value.name = "prob"
  )
  a = df[-ncol(df)]
  a = apply(a, 2,as.character)
  a = apply(a, 2,as.numeric)

  # anal is rank by probability
  anal = sapply(c(1:nrow(a)), function(x) rank(-a[x,],ties.method = "min"))
  anal = t(anal)

  Summary_t = anal %>%
    reshape2::melt()%>%
    group_by(Var2,value) %>%
    dplyr::summarize(n = n())


  #####important revision is needed -> DONE
  Rey = dcast(data = Summary_t,
             formula  = Var2~value,
             fun.aggregate = sum,
             value.var = "n"
)
  vec = c("Var2", as.character(1:length(unique(V_list$GENOTYPE))))

  #???????????????? No idea what is this doing
  if(length(setdiff(vec,  names(Rey)))>0){
    vecc = setdiff(vec,  names(Rey))
    for(p in 1:length(vecc)){

      Rey$new = 0
      names(Rey)[which(names(Rey)=="new")] = vecc[p]
    }
  }
  Rey = subset(Rey, select = vec)


  return(list(mat = df1, rankk = Rey))

}
