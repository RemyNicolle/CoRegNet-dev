\name{masterRegulator}
\alias{masterRegulator,coregnet-method}
\alias{masterRegulator}

\title{
Identify phenotype related Master Regulators.
}
\description{
This function implements methods to identify Master Regulators based on a regulatory network and on various type of input representing the implication of sets of gene in a phenotype of interest. These are derived from Celine Lefebvre's algorithm MARINa (Lefebvre, 2010).

}
\usage{
masterRegulator(coregnet,targetGenes,method=c("set.overlap","merge.pvalues","list.enriched"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coregnet}{
A large scale co-regulatory network of type coregnet

}
  \item{targetGenes}{
A liste of target genes. Can be either given as a character vector containing the genes or a named numeric vector containing weights (\emph{e.g.} fold change) or p-values with gene as names.
}
  \item{method}{The method to use to find the Master Regulators of the input target genes. The default is \code{set.overlap}. See details.
  }
}
\details{
Three types of input can be used depending on the objective.

To identify Master Regulators of a given set of genes (from a pathway, a liste of differentially expressed genes etc ...), \code{MasterRegulatorInference} simply needs a CoRegNet network object and a character vector describing the Target Genes of interest. Fisher's exact test will be used to identify TF with the set of target genes that is the most specific to these genes of interest.

To identify Master Regulators of a phenotype of interest, the p-values of the comparison with a reference phenotype (using a moderate or unmoderate t-test) or the Fold change can be used with a combined Fisher's test or a Kolmogorov-Smirnov test to identify significant TF of these genes. 

}

\value{
The sorted list of TF in the input network with it's associated p-value.
}
\references{
Lefebvre C, Rajbhandari P, Alvarez MJ, Bandaru P, Lim WK, Sato M, Wang K, Sumazin P, Kustagi M, Bisikirska BC, Basso K, Beltrao P, Krogan N, Gautier J, Dalla-Favera R and Califano A (2010) A human B-cell interactome identifies MYB and FOXM1 as master regulators of proliferation in germinal centers. Molecular Systems Biology 6: 1-10
}
\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}

\examples{
# Dummy expression data and network
gexp=matrix(rnorm(2600,sd=3),ncol=100)
gexp=rbind(gexp,do.call(rbind,lapply(1:26,function(i){
tf = sample(1:26,4)
  return((gexp[tf[1],]+gexp[tf[2],] -gexp[tf[3],]-gexp[tf[4],] +rnorm(100,sd=3))/2)})))
dimnames(gexp)=list(c(letters,LETTERS),paste("s",1:100,sep=""))
GRN=hLICORN(gexp,TFlist = letters)


MR=masterRegulator(GRN,LETTERS[1:10])
head(MR)

exampleWeight = rnorm(26)
names(exampleWeight) = LETTERS[1:26]
MR=masterRegulator(GRN,exampleWeight,"list")
head(MR)

examplePvalue = 10^-(0.1+runif(26))
names(examplePvalue) = LETTERS[1:26]
MR=masterRegulator(GRN,examplePvalue,"merg")
head(MR)

}
