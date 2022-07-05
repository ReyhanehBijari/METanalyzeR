#' Sample dataset for the genotype ranking
#'
#' A dataset containing the genotype and the locations and phenotype.
#'
#' @format A dataframe with 1688 observations and 10 variables:
#' \describe{
#'    \item{X}{X}
#'    \item{GENOTYPE}{genotypes planted in a trial}
#'    \item{ENVIRONMENT}{location-year the crop genotype is planted in}
#'    \item{REPNO}{replications of each genotype in each environment}
#'    \item{PT}{phenotypic of phenotype}
#'    \item{YEAR}{year}
#'    \item{CHECK}{whether the genotype is one of the baselines or not}
#'    \item{LOCATION}{location the crop genotype is planted in}

#' }
"sample_data"

oat_data = metan::data_ge[,1:4]

names(oat_data)[which(names(oat_data)=="ENV")]  = "ENVIRONMENT"
names(oat_data)[which(names(oat_data)=="REP")] = "REPNO"
names(oat_data)[which(names(oat_data)=="GY")] = "PT"
names(oat_data)[which(names(oat_data)=="GEN")] =  "GENOTYPE"


oat_data$YEAR = 2020
#oat_data$CHECK = F
oat_data$LOCATION = oat_data$ENVIRONMENT
