#' MET Data Preparation
#' Normalizes the PT
#' @param data the input data file which is a data frame
#' @param LocationCol the location column of the data
#' @param YearCol the year column of the data
#' @param GenotypeCol the genotype column of the data
#' @param PTCol the phenotype/desired output of the data
#' @import dplyr
#' @import reshape2
#' @export
#' @return A data frame with the normalized PT
#'
#'
MET_data_preparation = function(data,LocationCol,YearCol=NULL,GenotypeCol,PTCol){

  if (is.null(YearCol)){
    New_Environment = LocationCol
    print("If your dataset contains data in more than 1 year, please enter your year column as an input to the MET_data_preparation function")
  } else {
    New_Environment = paste(LocationCol,YearCol)
  }


  if (is.null(data)){
    stop("Please enter you dataset as an input to the MET_data_preparation function")
  }

  if (is.null(LocationCol)){
    stop("Please enter the location column of your dataset as an input to the MET_data_preparation function")
  }

  if (is.null(PTCol)){
    stop("Please enter the desired phenotype column of your dataset as an input to the MET_data_preparation function")
  }

  if (is.null(GenotypeCol)){
    stop("Please enter the genotype column of your dataset as an input to the MET_data_preparation function")
  }


  data$ENVIRONMENT  = New_Environment
  data$PT  = PTCol
  data$GENOTYPE = GenotypeCol


  P_Env = data %>%
    group_by(ENVIRONMENT) %>%
    summarize(P_Env = mean(PT, na.rm = T))

  data = subset(data, select = c(ENVIRONMENT,GENOTYPE,PT))

  data  = merge(data, P_Env,by = c("ENVIRONMENT"))

  data$Corrected_PT = data$PT - data$P_Env
  result  = data %>%
    group_by(GENOTYPE, ENVIRONMENT) %>%
    mutate(M_Corrected_PT = mean(Corrected_PT))
  #result = result %>% dplyr::select(-CHECK)

  return(result)
}
