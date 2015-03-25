\name{regulators}
\alias{regulators,coregnet-method}
\alias{activators,coregnet-method}
\alias{repressors,coregnet-method}
\alias{targets,coregnet-method}
\alias{regulators}
\alias{activators}
\alias{repressors}
\alias{targets}

\title{
Interogate a coregnet object.
}
\description{
Query the network for regulators of specific targets and targets of specific genes.
}
\usage{
regulators(object, target = NULL, type = c("single", "coregulators"))
activators(object,target,type=c("single","coregulators"))
repressors(object,target,type=c("single","coregulators"))
targets(object,regulator=NULL,type=c("regulating","activating","repressing"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The network in the form of a \code{coregnet} object to query.
}
  \item{target}{
  The target gene to query.
}
\item{regulator}{
  The regulator to query.
}
  \item{type}{
  The type of regulation to obtain. Differs depending on the function used.
}
}

\value{
For \code{regulators} if no target is given, returns integer vector with the number of targets for each regulators of the network. Given a non null vector,  a vector of the union of the regulators of all the genes is returned.
For \code{activators} and \code{repressors} the behavior is similar except that a target gene is needed.
If \code{type = "coregulators"} then only the regulators, activators or repressors which are found to be co-regulators, co-activators or co-repressors of the target genes are given.


\code{targets} with no given regulator returns a character vector of all the target genes in the network. Specifying a vector of regulators will return a vector of the union of the targets of all these regulators. The type of regulation can be specified to return only the activated or repressed targets.

}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}

\examples{
acts=apply(matrix(rep(letters[1:4],7),nrow=2),2,paste,collapse=" ")[1:13]
reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
co=coregnet(grn)

regulators(co)
regulators(co,"A")
regulators(co,"A","coregulators")


activators(co,"A")
activators(co,"A","coregulators")

repressors(co,"A")
repressors(co,"A","coregulators")

targets(co)
targets(co,"a")
targets(co,"a","reg")
targets(co,"a","act")
targets(co,"a","rep")
targets(co,c("a","b"),"act")


}
