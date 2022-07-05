#' RM Correction
#'
#' It removes the effect of the locations
#'
#' @param data the preprocessed file from the output of function
#' @export
#' @return Returns the data with a new feature that describes the corrected PT
#'
rm_correction = function(data){

  #data = data[which(data$EXPERIMENT==exp),]
  model = lm(PT~RM,data=data)

  sl = model$coefficients[2]

  int = model$coefficients[1]

  correction = mean(data$PT)-(int+data$RM*sl)

  data$CORRECTED_PT_RM = data$PT + correction

  return(data)
}
