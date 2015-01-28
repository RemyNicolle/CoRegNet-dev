\name{hLICORN}
\alias{hLICORN}

\title{
Hybrid Learning of co-operative regulation network.
}
\description{
Parallelized inference of co-regulatory network from gene expression only.
}
\usage{
hLICORN(numericalExpression, discreteExpression = 
discretizeExpressionData(numericalExpression), 
TFlist, GeneList = setdiff(rownames(numericalExpression), TFlist),
parallel = c("multicore", "no", "snow"), cluster = NULL,
minGeneSupport=0.1,minCoregSupport = 0.1,
maxCoreg=length(TFlist),
searchThresh=1/3,nGRN=100,
verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numericalExpression}{
A numerical Matrix containing the expression of genes and of
transcription factors that will be used to inferred the network.
Rownames should contain the Gene names/identifiers. 
Samples should be in columns but Colnames are not important.
The data will be gene-centered but not scaled.

}
  \item{discreteExpression}{
  
  optional. Should be in exactly the same format as \code{numericalExpression} (dimensions, colnames and rownames) and should contain value only in {-1,0,1} with -1 for under-expressed, 0 for no change and 1 for over expressed. For default value see details.
  

}
  \item{TFlist}{
A character vector containing the names of the genes that
are designated as Transcription Factor or Regulators. 
These should be contained in the \code{rownames} of 
the \code{numericalExpression} argument. Use \code{data(HumanTF)} for gene symbols of Human transcription factor
}
  \item{GeneList}{
optional. The list of genes for which Gene Regulatory Networks
should be inferred. Should be in the rownames of the expression data. 
If not provided will be taken as all the genes in the rownames
of the expression data that are not annotated as TF in \code{TFlist}.
}
  \item{parallel}{
optional. 
The type of parallel method to use to run the process on several threads.
Should be a value in "no" for single thread calculation, 
"multicore" to use several threads on the same local machine or "snow" 
to use the snow package in which case the cluster argument
should contain a cluster object (e.g. \code{makeCluster})
}
  \item{cluster}{
A object of type cluster from the snow package which describes 
a cluster of computer to run the inference of each gene in parallel.
Needed only when parallel is set to snow.
}
  \item{minGeneSupport}{
A float between 0 and 1. Minimum number of samples in which a gene has to have non-zero discretized expression value to be considered for regulatory network inference.
}
 \item{minCoregSupport}{
A float between 0 and 1. Minimum number of samples in which a set of co-regulators have the same non-zero discretized expression value to be considered as a potential set of co-regulator. Default is 0.1 . Can be increased in case of limitations on memmory.
}
 \item{maxCoreg}{
A integer. Maximum size of co-regulator to consider. Default is set to the number of TF. Can be decreased in case of limitations on memmory.
}
 \item{searchThresh}{
A float between 0 and 1. Minimum proportion of sample in which a gene and a set of co-regulators must have the same non-zero discretized value in order to be considered as a potential co-regulator of the gene.
}
 \item{nGRN}{
  if NA, takes only the best GRN models.
}
  \item{verbose}{
Sets wether information will appear in console during computation.
}
}

\details{
A parallelized implementation of the h-LICORN alogrithm.
The inference of the network can be run in parallel using
either a single machine or on a cluster using the make cluster in the
\code{parallel} package (now a base R package).

The two required inputs are the continuous gene expression matrix
and the list of transcription factor in the rows of the
dataset. For Human, a list of TF is included as default data set
using \code{data(HumanTF)} for symbols or \code{data(HumanTF_entrezgene)} 
for entrez gene IDs.

In case of memmory overflow, a higher minCoregSupport (0.2 to 0.4) will perform more shallow search of co-regulators but will limit the quantity of memmory used.

In some cases, in particular for very large datasets (several hundreds and thousands of samples), a lower minCoregSupport (down to 0.01) can be used for deep search.
}



\value{

A object of type CoRegNet with gene regulatory network (GRN)
containing several solution per gene and specifying the
co-regulation interactions, the same network in the form
of an adjacency list (adjacencyList) 
and the inferred co-regulators (coRegulators). 

}
\references{
Elati M, Neuvial P, Bolotin-Fukuhara M, Barillot E, Radvanyi F and Rouveirol C (2007) LICORN: learning cooperative regulation networks from gene expression data.
Bioinformatics 23: 2407-2414

Chebil I, Nicolle R, Santini G, Rouveirol C and Elati M (2014) Hybrid Method Inference for the Construction of Cooperative Regulatory Network in Human. IEEE Trans Nanobioscience
}
\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}
\note{

The regulatory network inference process is done in several steps.
The speed of the inference is dependant on several factors. The number of genes, the density of the discrete expression matrix (number of 1 and -1).
}
\seealso{
 \link{discretizeExpressionData},\link{coregnet-class}
}
\examples{
# Dummy expression data
gexp=matrix(rnorm(2600,sd=3),ncol=100)

gexp=rbind(gexp,do.call(rbind,lapply(1:26,function(i){
tf = sample(1:26,4)
  return((gexp[tf[1],]+gexp[tf[2],] -gexp[tf[3],]-gexp[tf[4],] +rnorm(100,sd=0.5))/2)})))
dimnames(gexp)=list(c(letters,LETTERS),paste("s",1:100,sep=""))

## Simple example of network inference
dummyNet=hLICORN(gexp,TFlist = letters)

## Infer a network only on a subset of genes
subgene = unique(dummyNet@GRN$Target)[1:2]
dummyNet=hLICORN(gexp,TFlist = letters,GeneList=subgene)

## Discretize data based on a set of reference samples (here 10 first)
discexp = discretizeExpressionData(gexp,refSamples=1:10)
dummyNet=hLICORN(gexp,TFlist = letters,discreteExpression=discexp)


## The network can be queried using the following functions
# returns the hub regulators
regulators(dummyNet)
# get the regulators of a given gene
regulators(dummyNet,"A")
activators(dummyNet,"A")
targets(dummyNet)
targets(dummyNet,"b")


# or transformed into a data.frame
coregnetToDataframe(dummyNet)



}
\keyword{ Network inference }
