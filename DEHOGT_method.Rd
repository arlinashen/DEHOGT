\name{Differentially Expressed Hheterogeneous Overdispersion Gene Test}
\alias{DEHOGT}
\title{DEHOGT Documentation
}
\description{
This function implements a generalized linear model (GLM) approach to identify differentially expressed (DE) genes between two treatment groups. It supports both quasi-Poisson and negative binomial distributions for read count data.

}
\usage{
glm_func <- function(data, treatment, covariates = NULL, norm_factors = NULL, dist = "qpois", padj = TRUE, pval_thre = 0.05, l2fc = FALSE, l2fc_thre = NULL, num_cores = 1)
}
\arguments{
 \item{data}{A numeric matrix of N rows (genes) by S columns (samples) containing read count data.}
  \item{treatment}{An integer vector of length S indicating the treatment group for each sample.}
   \item{covariates}{An optional matrix of gene-wise covariates. Default is NULL.}
  \item{norm_factors}{(Optional) A numeric vector of normalization factors for each sample. Defaults to NULL, indicating that no normalization is applied.}
  \item{dist}{Character string specifying the distribution to use. Can be "qpois" for quasi-Poisson or "negbin" for negative binomial. Defaults to "qpois".}
  \item{padj}{Logical indicating whether to adjust p-values for multiple testing using the Benjamini-Hochberg procedure. Defaults to TRUE.}
  \item{pval_thre}{Numeric specifying the p-value threshold for identifying DE genes. Defaults to 0.05.}
  \item{12fc}{Logical indicating whether to filter results based on log2 fold change. Defaults to FALSE.}
  \item{12fc_thre}{Numeric specifying the log2 fold change threshold for identifying DE genes if l2fc is TRUE. Defaults to NULL.}
  \item{num_cores}{The number of CPU cores to use for parallel computing. Default is 1.}

}
\details{
This generalized linear model (GLM) approach is for identifying differentially expressed (DE) genes across two treatment groups. The function is adept at handling both quasi-Poisson and negative binomial distributions for the analysis of read count data, making it versatile for various types of gene expression datasets. The function outputs a list containing indices of DE genes, their p-values, log2 fold changes, and the normalization factors used in the analysis.
}
\value{
\describe{
\item{DE_idx}{A binary vector indicating whether each gene is differentially expressed (1) or not (0) based on the specified criteria.}
\item{pvals}{A numeric vector of adjusted (or unadjusted) p-values for each gene.}
\item{log2fc}{A numeric vector of log2 fold changes for each gene.}

}


\author{
Yubai Yuan, Qi Xu, Arlina Shen, Annie Qu
}
\examples{
# Generate example data
data <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
treatment <- sample(0:1, 10, replace = TRUE)
covariates <- matrix(rnorm(100), nrow = 100, ncol = 5)
# Run glm_func with parallel computing using 4 cores
result <- glm_func(data, treatment, covariates, num_cores = 4)


}
\references{
Yuan Y, Xu Q, Wani A, Dahrendor J, Wang C, Donglasan J, Burgan S, Graham Z, Uddin M, Wildman D, Qu A. Differentially Expressed Heterogeneous Overdispersion Genes Testing for Count Data. bioRxiv [Preprint]. 2023 Feb 22:2023.02.21.529455. doi: 10.1101/2023.02.21.529455. PMID: 36865247; PMCID: PMC9980115.}
