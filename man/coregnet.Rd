\name{coregnet}
\alias{coregnet}

\title{
Initialize a co-regulatory network object.
}
\description{
Given a previously built network encoded in a simple data.frame, initialize a co-regulatory network object \code{coregnet}. This can be used on networks inferred from gene expression data by other algorithms (ARACNe, GENIE3 ...) but also on ChIP-seq data or on network based on Transcription Factor Binding Site analysis. If the input network does not contain a dscription of the type of regulation (activation, inhibition) an expression data set is needed.
}
\usage{
coregnet(GRN, expressionDATA = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GRN}{
  A data.frame containing the network. A two column data.frame should contain the regulator in column one and the target gene in column two. A three column data.frame should contain the target in column one, all the co-activators and co-repressors seperated by a space in column two and three respectively.
}
  \item{expressionDATA}{
If the input GRN is only a two column data.frame, the expression data is needed to set each regulator of a target gene as an activator or a repressor.
}
}
\details{
In the case of a two column data.frame, the pearson correlation coefficient is used to determine wether a given regulator is an activator (R-squared >= 0) or a repressor (R-squared < 0).

In order to import a network and initialize a new coregnet object the input data.frame should have three columns, the first containing the target genes, then the activators and finally the inhibitors. Target genes can be present in multiple lines. Several regulators can be present in column 2 and 3 if they are separated by a space. Below is an example of a toy network :

A;R1 R2;R3 R4

A;R1 R5;R3 R6

B;R1 R2 R3;NA

...
}
\value{
Returns a \code{coregnet} object.
}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}
\seealso{
\code{\link{hLICORN}}
}
\examples{
acts=apply(matrix(rep(letters[1:4],4),nrow=2),2,paste,collapse=" ")
reps=apply(matrix(rep(letters[5:8],4),nrow=2),2,paste,collapse=" ")
grn=data.frame("Target"= LETTERS[1:16] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(16),stringsAsFactors=FALSE)
co=coregnet(grn)

}

