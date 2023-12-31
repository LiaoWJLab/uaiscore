% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_cox_model.R
\name{predict_cox_model}
\alias{predict_cox_model}
\title{Function to predict risk scores and feature importance using a random forest (RF) model in a Cox regression setting}
\usage{
predict_cox_model(
  sur_model,
  eset_new,
  pdata_new = NULL,
  feas = NULL,
  feas_genes = NULL,
  prefix = c("-", "\\\\:"),
  id_pdata_new = "ID",
  na_replace = 1
)
}
\arguments{
\item{sur_model}{Survival model object}

\item{eset_new}{New expression matrix data}

\item{pdata_new}{New patient metadata (default: NULL)}

\item{feas}{List of features (default: NULL)}

\item{feas_genes}{Genes associated with the features (default: NULL)}

\item{prefix}{Prefix for column names (default: c("-", "\\:"))}

\item{id_pdata_new}{Identifier for patient metadata (default: "ID")}

\item{na_replace}{NA replacement value (default: 1)}
}
\description{
Function to predict risk scores and feature importance using a random forest (RF) model in a Cox regression setting
}
\author{
Dongqiang Zeng
}
