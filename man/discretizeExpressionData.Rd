\name{discretizeExpressionData}
\alias{discretizeExpressionData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Three-value discretization of gene expression data.
}
\description{
Pre-process step to transform log2 numerical expression data into 
a three value categorical : over-expression (+1), under-expression (-1)
and no change (0).
}
\usage{
discretizeExpressionData(numericalExpression, threshold = NULL,
refSamples = NULL, standardDeviationThreshold = 1)

discretizeExpressionData(numericalExpression, threshold = NULL,
refSamples = NULL, standardDeviationThreshold = 1)

}

\arguments{
  \item{numericalExpression}{
A matrix of continuous log2 gene expression data with genes
in rows and samples in column.
}
  \item{threshold}{
A numeric value used as a fixed fold change threshold for brut discretization. Can
be NULL if a standardDeviationThreshold is not.
}
  \item{refSamples}{
A vector of column names used as a set of reference samples to be used to compute fold changes.
Can be NULL if the input data is already centered on a reference
values (tested by the presence of negative values) or if the mean
if each gene should be used to center (not scale) each expression values.
}
  \item{standardDeviationThreshold}{
The multplicator of the whole data set standard deviation to
be used as the threshold.
}
}
\details{
Given a continuous log2 gene epression matrix this function aims
at producing a matrix of discretized expression values.
The numerical data must be in some form of fold change to compare
each value of a gene in a sample with a reference value of the same gene.
The behavior of the function will therefore depend on the form of the input data.

Given a matrix with negative values, the function will 
consider that the data is already in the right format and
will simply apply a hard threshold to discetize the data.

Given a matrix with only positive values, which is the case
for normalized RNAseq or single color microarrays, the function
will center each gene based on it's mean expression in all samples
or based on the mean of expression of a set of reference sample
(normal samples in a study of a particular disease for example).

In either case, the threshold will be used to transform the data
in +1s if the value of a gene in a sample is above or equal to the
threshold, -1s if the value is below the the negative value of the threshold
and 0 otherwise.

The default is to compute a threshold based on the overall distribution of the
numerical values in the dataset. This was choosen over a default fold change example
(usually 1 or 2 corresponding to a two-fold or four-fold increase/decrease)
after observing a large difference between technologies. However, the choice
between a simple hard fold change threshold or a threshold as a multiplicator 
of the global standard deviation remains.


}
\value{
A matrix of integers with the same number of rows and the same number of column as the input numericalExpression. Values in the output matrix are in {-1,0,1}.
The reference samples are removed if given.
}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}
\examples{
# Use mean of each gene as a reference
expression=matrix(2*rnorm(200),nrow=2,dimnames=list(paste("gene",1:2,sep=""),paste("sample",1:100,sep="")))
discExp=discretizeExpressionData(expression)
boxplot(expression~discExp,xlab="Discrete values",ylab="Continuous values")
pie(table(discExp))
discExp=discretizeExpressionData(expression,standardDeviationThreshold=2)
pie(table(discExp))
discExp=discretizeExpressionData(expression,threshold=1)
pie(table(discExp))

# Use of reference sample
expression=matrix(2*rnorm(200),nrow=2,dimnames=list(paste("gene",1:2,sep=""),paste("sample",1:100,sep="")))
discExp=discretizeExpressionData(expression,refSamples=1:10)
boxplot(expression~discExp,xlab="Discrete values",ylab="Continuous values")
pie(table(discExp))
discExp=discretizeExpressionData(expression,standardDeviationThreshold=2,refSamples=1:10)
pie(table(discExp))
discExp=discretizeExpressionData(expression,threshold=1,refSamples=paste("sample",1:10,sep=""))
pie(table(discExp))

}

