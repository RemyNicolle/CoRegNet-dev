\name{addEvidences}
\alias{addEvidences,coregnet-method}
\alias{addCooperativeEvidences,coregnet-method}
\alias{addEvidences}
\alias{addCooperativeEvidences}

\title{
Integration of regulatory evidences.
}
\description{
These functions can be used to integrate external data on gene regulation and co-regulation to enrich an inferred regulatory network.
Additional regulatory data sets can include : ChIP-seq data, ChIP on chip, target gene promoter analysis for Transcription Factor Binding Site models, epigenetic marks or even networks inferred by other methods. 
Cooperative regulation data sets can include : Protein interactions, significant binding site overlap, co-localization.
These additional evidence are added to enrich the network and the enrichment of the inferred interaction is tested. 
}
\usage{
addEvidences(object,...)
addCooperativeEvidences(object,...)
}

\arguments{
  \item{object}{
    A regulatory network inferred by the hLICORN function or with an externally inferred network importing through the CoRegNet class.
}
    \item{...}{
    A single or several data sets of regulatory data as pairs of Transcription Factor to Target Gene regulation interactions or of pairs of transcription factors for cooperative evidences. These should be given in the form of a two column data.frame or matrix. The characters in the input dataset should correspond to the characters of the regulators and the genes in the network (accessible through \link{regulators} \link{targets} respectively).}
}

\details{

A single or several datasets of regulatory interactions can be added and enrich an inferred regulatory network. Below is an example of the format of the input data
  A character matrix or data.frame with two columns the first for target genes and the second for Tanscription Factor.
  
     [,1] [,2]
     
[1,] "TF1"  "Gene1" 

[2,] "TF1"  "Gene2" 

[3,] "TF2"  "Gene1" 
  
  ...
  
}




\value{
A Gene Regulatory Network in a CoRegNet objects with additional informations in the core object.
}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}
\note{
The variables used as the input are kept and will be used later on to refine the network.
The names of target genes and of TF need to be identical to the ones in the expression dataset.
}

\seealso{
 \link{refine}
 \link{hLICORN}
}
\examples{

#Creating a synthetic network. Upper case letters are target genes, lower case letters are regulators
acts=apply(matrix(rep(letters[1:4],4),nrow=2),2,paste,collapse=" ")
reps=apply(matrix(rep(letters[5:8],4),nrow=2),2,paste,collapse=" ")
grn=data.frame("Target"= LETTERS[1:16] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(16),stringsAsFactors=FALSE)
co=coregnet(grn)

#Creating dummy regulatory evidence data in a data.frame
evidence1=unique(data.frame(tf=sample(letters[1:8],50, replace=TRUE),target=sample(LETTERS[1:16],50, replace=TRUE),stringsAsFactors = FALSE))

evidence2=unique(data.frame(tf=sample(letters[1:8],50, replace=TRUE),target=sample(LETTERS[1:16],50, replace=TRUE),stringsAsFactors = FALSE))


#Creating dummy cooperative vidence data in a data.frame
coregevidence1=unique(data.frame(tf1=sample(letters[1:8],50, replace=TRUE),target=sample(letters[1:8],50, replace=TRUE),stringsAsFactors = FALSE))



enrichco=addEvidences(co,evidence1,evidence2)
enrichco=addCooperativeEvidences(enrichco,coregevidence1)


}

