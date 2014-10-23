\name{display}
\alias{display}

\title{
Display a shiny interactive web interface to 
}

\description{
Launches a shiny webpage for interactive viewing and analysis of the co-regulation network using a javascript cytoscape network. 
}

\usage{
display(coregnetwork, expressionData = NULL, TFA = NULL, alterationData = NULL, clinicalData = NULL, TFnotes = NULL, allTFplot = .heatplot, oneTFplot = .tfPlot)
}

\arguments{
  \item{coregnetwork}{
  A coregnet object
}
  \item{expressionData}{
  A matrix or data.frame object with named columns (samples) and named rows (genes) containing gene expression data
}
  \item{TFA}{
  A matrix or data.frame object with named columns (samples) and named rows (genes) containing transcription factor activity data such as the data obtained using the \code{regulatorInfluence} function. (unnecessary if expression data is given, makes the function run faster)
}
  \item{alterationData}{
  optional. A matrix or data.frame object with named columns (samples) and named rows (genes) containing gene alteration data
}
  \item{clinicalData}{
  optional. Either a list or a factor describing clinical information about samples. A list must be named and each entry should contain a vector of samples. 
}
  \item{TFnotes}{
 optional. A factor describing TFs.
}
  \item{allTFplot}{
  A function ploting information by default. Default functions are implemented.
}
  \item{oneTFplot}{
  A function ploting information about a TF which will be used when a single node is selected on the network. Default function are implemented.
}
}



\value{
Does not return anything.
}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}



\examples{
acts=apply(matrix(rep(letters[1:4],7),nrow=2),2,paste,collapse=" ")[1:13]
reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
co=coregnet(grn)
samples= paste("S",1:100,sep="")
expression=matrix(rnorm(3400),ncol=100)
dimnames(expression) = list(c(grn$Target,names(regulators(co))),samples)

TFA = regulatorInfluence(co,expression,minTarg=4)
colnames(TFA) = samples
if(interactive()){
 display(co,TFA=TFA,expressionData=expression)
}

CNA =matrix( sample(-2:2,800,replace=TRUE),ncol=100)
dimnames(CNA) = list(names(regulators(co)),samples)
if(interactive()){
  display(co,TFA=TFA,expressionData=expression,alteration=CNA)
}

clinicGrp = factor(paste("grp",sample(1:3,100,replace=TRUE),sep=""))
names(clinicGrp) =samples
if(interactive()){
 display(co,TFA=TFA,expressionData=expression,alteration=CNA,clinicalData=clinicGrp)
}
}

\keyword{shiny}
\keyword{cytoscape}
\keyword{interactive}
\keyword{network}