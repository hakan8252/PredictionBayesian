% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PredictionBayesian.R
\name{Cross_Validation}
\alias{Cross_Validation}
\title{Cross Validation Method for Bayesian Network Models}
\usage{
Cross_Validation(
  fold,
  data.tt,
  target,
  input,
  num.iter,
  str.algorithms = NULL,
  structure.list = NULL
)
}
\arguments{
\item{fold}{an integer value shows the fold number of cross validation.}

\item{data.tt}{a data frame containing the variable(s) in the Bayesian network model(s). Values of the variables
in the data frame are assumed to be discrete.}

\item{target}{target variable(s) whose actual values and posterior distributions will be obtained.}

\item{input}{input variable(s) to be used when calculating the posterior distribution of target variables.}

\item{num.iter}{an integer that shows number of query variables which will be used during calculation of
posterior distribution of target variable. Default value is `1`.}

\item{str.algorithms}{a list containing different algorithms \code{\link[bnlearn]{structure-learning}} that are used to learn
different Bayesian networks. The default value is `NULL`.}

\item{structure.list}{If there is pre-defined structure of Bayesian network, parameters of it are learned. The default value is `NULL`.}
}
\value{
a list containing AUC and ROC curve results. List has two index. First index consists of data frame containing AUC results of
        all folds. Second index consists of ROC curve results of all folds.
}
\description{
This function is used to perform k-fold cross validation method for a Bayesian network models.
}
\details{
by using \code{\link{Predictions}} and \code{\link{Actual.Values}} functions posterior distribution and actual values of target
         variable(s) obtained. With the \code{\link{Model_Learning}} function Bayesian network model(s) is learned. In the last step
         ROC curves and AUC values are obtained through \code{\link{ROC_Calculation}} and model performances are evaluated. K-fold
         cross validation method is applied.
}
\references{
M. Stone, “Cross-Validatory Choice and Assessment of Statistical Predictions,” J. R. Stat. Soc. Ser. B, vol. 36, no. 2, pp. 111–133, 1974.
}
\seealso{
\code{\link{information.algorithm}} \code{\link[bnlearn]{bn.cv}}
}
