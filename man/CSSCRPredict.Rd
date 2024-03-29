% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSSCRPredict.R
\name{CSSCRPredict}
\alias{CSSCRPredict}
\title{Make prediction with the CSCR model}
\usage{
CSSCRPredict(
  old_data_x = NULL,
  new_data_x = NULL,
  old_loading_x,
  old_loading_y,
  mean_x,
  new_assign = NULL
)
}
\arguments{
\item{old_data_x}{a matrix indicating the training predictor set}

\item{new_data_x}{a vector indicating the test predictor set of a new observation}

\item{old_loading_x}{a matrix indicating the loading matrix associated with the training data set}

\item{old_loading_y}{a vector indicating the regression coefficients associated with the training data set}

\item{mean_x}{the proportion of cluster differences explained by the mean-level differences}

\item{new_assign}{whether the cluster assignment of the new observation is known}
}
\value{
two items indicate the prediction of the outcome and of the cluster assignment for the new entry
}
\description{
This function allows modellers to predict an new observations based on the information with the existing CSSCR analyis
}
