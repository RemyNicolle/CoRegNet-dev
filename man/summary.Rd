\name{summary}
\alias{summary,coregnet-method}
\alias{dim,coregnet-method}
\alias{length,coregnet-method}
\alias{show,coregnet-method}
\alias{print,coregnet-method}
\alias{coregnetToDataframe,coregnet-method}
\alias{coregnetToList,coregnet-method}
\alias{coregnetToList.CoRegNet}
\alias{coregnetToList.coregnet}
\alias{summary}
\alias{coregnetToList}
\alias{coregnetToDataframe}
%\alias{coerce,coregnet-method}
%\alias{coerce,coregnet-method,data.frame}
%\alias{coerce,coregnet-method,list}
%\alias{as.data.frame,coregnet-method}
%\alias{as.list,coregnet-method}

\title{Summaries and info coregnet}
\description{
Several functions to print and view info about the network enclosed in a coregnet object.
}
\usage{
summary(object,...)
\S4method{show}{coregnet}(object)
\S4method{dim}{coregnet}(x)
\S4method{length}{coregnet}(x)
\S4method{print}{coregnet}(x)
\S4method{coregnetToDataframe}{coregnet}(network)
\S4method{coregnetToList}{coregnet}(network)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{network}{
a coregnet network object.
}
  \item{object}{
a coregnet network object.
}
\item{x}{
a coregnet network object.
}

  \item{...}{
unused argument
}
}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}


\examples{
regs=sample(letters,7)
grn=data.frame("Target"= LETTERS ,"activators"= sample(rep(regs,4))[1:26],"repressors"= sample(rep(regs,4))[1:26],stringsAsFactors=FALSE)
co=coregnet(grn)
print(co)
length(co)
dim(co)
co
coregnetToDataframe(co)
coregnetToList(co)
}

