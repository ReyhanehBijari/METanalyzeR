% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MET_data_preperation.R
\name{MET_data_preparation}
\alias{MET_data_preparation}
\title{MET Data Preparation
Normalizes the PT}
\usage{
MET_data_preparation(data, LocationCol, YearCol = NULL, GenotypeCol, PTCol)
}
\arguments{
\item{data}{the input data file which is a data frame}

\item{LocationCol}{the location column of the data}

\item{YearCol}{the year column of the data}

\item{GenotypeCol}{the genotype column of the data}

\item{PTCol}{the phenotype/desired output of the data}
}
\value{
A data frame with the normalized PT
}
\description{
MET Data Preparation
Normalizes the PT
}
