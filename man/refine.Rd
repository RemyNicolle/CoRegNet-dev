\name{refine}
\alias{refine,coregnet-method}
\alias{refine}

\title{
    Refine an inferred regulatory network using external evidence.
}
\description{
  Refines the inferred network using the integrated external evidences added by \link{addEvidences} or \link{addCooperativeEvidences}. Several strategies can be applied depending on the number and type of added evidences. These include supervised and unsupervised processes to use all the integrated data set with the inferred network and select the best Gene Regulatory Networks (GRN). Can also be used when no additional dataset has been integrated.
}
\usage{
refine(object,GRNselection=c("best","maximize","threshold"),
  integration=c("unsupervised","supervised"),
  referenceEvidence=NULL,evidenceToMaximize="R2",threshold=NULL,verbose=TRUE)
}

\arguments{
  \item{object}{
   A regulatory network inferred by the hLICORN function which can also have been enriched by external regulatory data sets using   \link{addEvidences} or \link{addCooperativeEvidences}.
}
  \item{GRNselection}{
  The type of Gene Regulatory Network (GRN) selection method to apply. Default is to select the best regulatory model per gene. See details for other possibilities.
}
  \item{integration}{
Defines the method to merge all the available data sets into a single score which will define the quality of a given GRN. Default is unsupervised. See details.
}
  \item{referenceEvidence}{
To be specified when using the supervised integration method. Specifies one of the integrated data set as a Gold standard to learn the best weight of each evidences to maximize the number of reference evidence in the final network.
}

  \item{evidenceToMaximize}{
  To be specified when using the \code{maximize} GRN selection scheme. Instead of selecting one GRN per gene, this method will choose a threshold for the merged score that will maximize the number of interaction from a given evidence data set. When using the \code{supervised} integration method, the default is to use the \code{referenceEvidence} to be maximized.
}
  \item{threshold}{
  When the automatically choosen threshold is not satisfactory, a user given threshold can be applied.
  }
  
  \item{verbose}{
  If set to TRUE (the default) sends messages at each step of the process.
  }
}

\details{
This function implements several strategies to select the best large scale regulatory network. Depending on the number and type of added evidences the strategies and some recommandations are detailed below.

The first step of the refinement is the integration of the different external evidences into a \emph{merged score}. If no evidence data set has been added, the score given by the inference algorithm is used by its own (an adjusted R2). In the \code{unsupervised} method, the default, the merged score is simply the mean of each of the evidences, including the score given by hLICORN. For the \code{supervised} process, a weight is given to each of the evidences which is learned using a generalized linear model which will be fitted to predict the regulatory interactions of a user defined reference evidence data set.

These two (\code{unsupervised} or \code{supervised}) integration methods are derived from the network learning process used by the modENCODE consortium (Marbach et al, 2012). Once the merged score is obtained, the default is to select the best GRN per gene. However, another possibility is to select all the "good" networks by choosing a threshold on the merged score. This can be done either by a user defined threshold between 0 and 1 or by choosing automatically a threshold that will maximize the interactions originating from a user defined evidence data set.

The default behavior of the function is to integrate the data set in an unsupervised way and selecting the best GRN per gene.

It is recommanded that when no additional data has been integrated and the selection of the network only needs to be done based on the inference score, a bootstraped regression coefficient, then the simplest strategies is to use the default parameters which will select the best GRN per gene.

}

\value{
A coRegNet object specifying a refined large scale co-regulatory network.
}



\references{
Marbach D, Roy S, Ay F, Meyer PE, Candeias R, Kahveci T, Bristow CA & Kellis M (2012) Predictive regulatory models in Drosophila melanogaster by integrative inference of transcriptional networks. Genome Research 22: 1334-1349
}
\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}

\seealso{
\link{addEvidences} and \link{addCooperativeEvidences}
}
\examples{
#Dummy network and evidence data examples
acts=apply(rbind(rep("z",14),matrix(rep(letters[1:4],7),nrow=2)),2,paste,collapse=" ")[1:13]
reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
GRN=coregnet(grn)

tfs=letters
genes=LETTERS
evidence1=unique(data.frame(tf=sample(tfs,100,replace=TRUE),target=sample(genes,100,replace=TRUE),stringsAsFactors =FALSE))
evidence2=unique(data.frame(tf=sample(tfs,100,replace=TRUE),target=sample(genes,100,replace=TRUE),stringsAsFactors =FALSE))
evidence3=unique(data.frame(tf=sample(tfs,100,replace=TRUE),target=sample(genes,100,replace=TRUE),stringsAsFactors =FALSE))

GRNenrich=addEvidences(GRN,evidence1,evidence2,evidence3)
print(GRNenrich)



unsupervisedNet=refine(GRNenrich)
supervisedNet=refine(GRNenrich,integration="sup",referenceEvidence="evidence1")

# The following usually gives poor results...
#supervisedNet=refine(GRNenrich,integration="sup",referenceEvidence="evidence1",evidenceToMaximize="evidence1",GRNselection="maximize")


}
