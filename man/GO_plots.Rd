\name{GO_plots}
\alias{GO_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform GO enrichment for a set of genes, and create plot
%%  ~~function to do ... ~~
}
\description{
This function takes as argument a vector of gene ENTREZ IDs that are significantly differentially expressed between two conditions, runs a GO enrichment analysis (MF, BP & CC) and produces ggplot2 objects for visualization of the results. In the ggplot the top 30 (or if less than 30, all) GO terms (by p.adjust value) are displayed.
}
\usage{
GO_plots(sigGenes, name="name")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigGenes}{
  A character vector containing gene ENTREZ IDs.
}
  \item{name}{
  A character string to be added to the ggplot2 objects.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
