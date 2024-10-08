% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_heatmap.R
\name{ig_scaleV}
\alias{ig_scaleV}
\title{return a vector of scaled sizes for a vector of verticies}
\usage{
ig_scaleV(vals, sf = 5, minVal = 2)
}
\arguments{
\item{vals}{values associated with verticies.  Can be number of sub-net nodes (members) or degree (number of edges))}

\item{sf}{scaling factor, so sf=5 means that the maximum vertix will have cex=5)}

\item{minVal}{minVal}
}
\value{
igraph object
}
\description{
return a vector of scaled sizes for a vector of verticies
}
