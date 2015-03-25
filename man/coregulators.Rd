\name{coregulators}
\alias{coregulators,coregnet-method}
\alias{coregulators}

\title{
Extract all co-regulators of a regulatory network.
}
\description{
Based on the frequency and specificity of co-regulation, this functions extracts from a coregnet network all the cooperative regulators.
}
\usage{
coregulators(object, maxcoreg = 2, verbose = TRUE, minCommonGenes = ifelse(maxcoreg == 2, 1, 10), adjustMethod = "fdr", alpha = 0.01)
}

\arguments{
  \item{object}{
 The coregnet object containing the co-regulatory network
}
  \item{maxcoreg}{
  The maximum size of co-regulator sets to extract. The default is 2 to extract pairs of co-regulators
}
  \item{verbose}{
To output information during the extraction process. By default does not ouput anything except warnings.
}
  \item{minCommonGenes}{
The minimum number of commmon genes, co-regulated genes, to consider two regulators for further co-regulation analysis. Default is 1 for pairs and 10 for bigger sets.
}
  \item{adjustMethod}{
The p-value adjustment method to extract significant pairs of co-regulators. Default is FDR correction. Anything that is accepted by the \code{\link{p.adjust}} method is fine.
}
  \item{alpha}{
The threshold to consider a pair of co-regulator significant (after pvalue correction). Default is 1\%
}
}



\value{
For \code{maxcoreg} set to 2, a data.frame with pairs of regulators in the two first columns (Reg1 and Reg2), the Support (the portion of coregulated gene networks in which this pair is found), nGRN for the number of times the pair was found as cooperative regulators in the net, a p-value of a FisherTest testing the specificty of the shared targets and the adjusted pvalue of this test.

For \code{maxcoreg} higher than 2, a two column data.frame with the set of co-regulators (collapsed in one character seperated by a space) and the Support.


}
\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}


\seealso{
\code{\link{regulators}},\code{\link{hLICORN}}
}

\examples{
acts=apply(rbind(rep("z",14),matrix(rep(letters[1:4],7),nrow=2)),2,paste,collapse=" ")[1:13]
reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
co=coregnet(grn)

coregulators(co)
coregulators(co,maxcoreg=3,minCommonGenes=3)

coregulators(co,adjustMethod="bonferroni")

coregulators(co,alpha=1)

}

