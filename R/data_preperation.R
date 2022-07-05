#' Data Preparation
#' Normalizes the PT
#' @param data the input data file which is a data frame
#' @import dplyr
#' @import reshape2
#' @export
#' @return A data frame with the normalized PT
#'
#'
data_preparation = function(data){

  #removing observation without harvest date
  #data = data%>%dplyr::filter(!HARVEST_DATE == "")
  if(!"ENVIRONMENT "%in% names(data)){
    data$ENVIRONMENT = paste(data$LOCATION,data$YEAR)
    #data$Breeder = substring(data$EXPERIMENT,5,6)
  }
  #finding average phenotype in each environment
  P_Env = data %>%
    group_by(ENVIRONMENT) %>%
    summarize(P_Env = mean(PT, na.rm = T))

  data  = merge(data, P_Env,by = "ENVIRONMENT")

  data$Corrected_PT = data$PT - data$P_Env
  result  = data %>%
    group_by(GENOTYPE, ENVIRONMENT) %>%
    mutate(M_Corrected_PT = mean(Corrected_PT))
    #result = result %>% dplyr::select(-CHECK)

  return(result)
}
