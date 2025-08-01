% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{neuronalSignature}
\alias{neuronalSignature}
\title{Compute ADAR activity with the neuronal signature}
\usage{
neuronalSignature(mat, object)
}
\arguments{
\item{mat}{Either a numeric expression matrix with genes as rows and samples
as columns or an ExpressionSet from a Seurat object, SingleCellExperiment or
SpatialExperiment}

\item{object}{Character string indicating the type of object, either bulk or
eset}
}
\value{
Matrix with ADAR activity scores
}
\description{
The function returns ADAR activity scores for each sample/spot/cell depending
on the starting expression matrix given as input. It employs a signature
derived from a dataset of hESC h9 cell line samples, which is ideal to
measure ADAR activity in a human neuronal context.
}
\examples{

# Load the regulon
data(neuronal_regulon)

# Load the example expression matrix
data(example_matrix)

# Compute ADAR activity
res <- neuronalSignature(example_matrix, object = "bulk")

}
