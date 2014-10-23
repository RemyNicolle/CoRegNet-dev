\name{regulatorInfluence}
\alias{regulatorInfluence,coregnet-method}
\alias{regulatorInfluence}

\title{
Regulator Influence, estimating the sample specific activity of Transcription Factors.
}
\description{
  Uses a network in the form of a coregnet object to compute regulatory influence to estimate the transcriptional activity of each regulators in each sample of the given expression data.
}
\usage{
  regulatorInfluence(object,expData,minTarg = 10,withEvidences=FALSE,addCoregulators=FALSE, is.scaled=FALSE) 
  }

\arguments{
  \item{object}{
  A network in the form of a coregnet object.
}
  \item{expData}{
  An expression data matrix or data.frame.
}
\item{addCoregulators}{
  Compute influence for coregulators with sufficient number of targets. Default to FALSE.
}
  \item{minTarg}{
  The minimum number of targets for a regulator to be considered for actvity prediction. Default set to 10. 
}
\item{withEvidences}{
  Use only the target genes which are validated by an external validation dataset (ChIP-seq for example). This is only possible if external evidence was added using \link{addEvidences}. Default set to False. 
}
\item{is.scaled}{
	Wether the input expression data is scaled, if not it will be.
  }

}

\value{
An N by R matrix with N columns the number of sample in the original expression data and R rows the number of regulators with sufficient targets to compute their influence.

The expression data is centered by default but not scaled.
}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}

\references{
Nicolle R, Elati M and Radvanyi F (2012) Network Transformation of Gene Expression for Feature Extraction. In pp 108-113.
}

\seealso{
\code{\link{hLICORN}} and \code{\link{coregnet-class}} to create the network.
\code{\link{addEvidences}} to add external evidences.
}
\examples{
acts=apply(matrix(rep(letters[1:4],7),nrow=2),2,paste,collapse=" ")[1:13]
reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
co=coregnet(grn)
samples= paste("S",1:100,sep="")
expression=matrix(rnorm(3400),ncol=100)
dimnames(expression) = list(c(grn$Target,names(regulators(co))),samples)

#Minimum number of targets is adjusted because of the small size of the network
TFA = regulatorInfluence(co,expression,minTarg=4)


}
