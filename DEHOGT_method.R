library(DESeq2)
library(foreach)
library(doParallel)
library(MASS)
library(DESeq2)
library(DESeqAnalysis)
library(ROSE)
library(edgeR)
library(multcomp)
library(MLmetrics)
library(VennDiagram)
library(pROC)
library(progress)
library(dplyr)
library(tidyverse)
library(ggplot2)
# Check again for BiocManager; if not present, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install DESeqAnalysis package from specified repositories
install.packages(
  pkgs = "DESeqAnalysis",
  repos = c(
    "r.acidgenomics.com",
    BiocManager::repositories()
  )
)
#' @param data A matrix of gene expression data where rows represent genes and columns represent samples.
#' @param treatment A vector specifying the treatment conditions for each sample.
#' @param covariates An optional matrix of gene-wise covariates. Default is NULL.
#' @param norm_factors An optional vector of normalization factors for each sample. Default is NULL, which assumes equal normalization factors.
#' @param dist The distribution family for the GLM. Can be "qpois" for quasi-Poisson or "negbin" for negative binomial. Default is "qpois".
#' @param padj Logical value indicating whether to adjust p-values using the Benjamini-Hochberg (BH) procedure. Default is TRUE.
#' @param pval_thre The threshold for identifying differentially expressed genes based on adjusted p-values. Default is 0.05.
#' @param l2fc Logical value indicating whether to consider log2 fold change for identifying differentially expressed genes. Default is FALSE.
#' @param l2fc_thre The threshold for log2 fold change in identifying differentially expressed genes. Default is NULL.
#' @param num_cores The number of CPU cores to use for parallel computing. Default is 1.
#'
#' @return A list containing:
#'   \item{DE_idx}{A logical vector indicating differentially expressed genes.}
#'   \item{pvals}{A numeric vector of p-values for each gene.}
#'   \item{log2fc}{A numeric vector of log2 fold changes for each gene.}
#'
#' @examples
#' data <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' treatment <- sample(0:1, 10, replace = TRUE)
#' covariates <- matrix(rnorm(100), nrow = 100, ncol = 5)
#' Run glm_func with parallel computing using 4 cores
#' result <- glm_func(data, treatment, covariates, num_cores = 4)

#' @import DESeq2
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @export
glm_func <- function(data, treatment, covariates = NULL, norm_factors = NULL, dist = "qpois", padj = TRUE, pval_thre = 0.05, l2fc = FALSE, l2fc_thre = NULL, num_cores = 1) {

  # Initialize parallel computing
  registerDoParallel(cores = num_cores)

  dds = DESeqDataSetFromMatrix(data, DataFrame(treatment), ~treatment)

  if (is.null(norm_factors)){
    norm_factors = rep(1, ncol(data))
  }

  result_GLM <- foreach(i = 1:nrow(data), .combine = rbind) %dopar% {

    genewise_data = data[i, ]
    genewise_dataframe = data.frame(read = genewise_data, treatment = treatment)

    if (!is.null(covariates)) {
      genewise_dataframe <- cbind(genewise_dataframe, covariates[i, ])
    }

    if (dist == "qpois"){
      model_formula <- as.formula(paste("read ~ treatment + offset(log(norm_factors))", if (!is.null(covariates)) paste("+", colnames(covariates), collapse = " + "), sep = " "))
      model = glm(formula = model_formula, family = quasipoisson(), data = genewise_dataframe)
    } else if (dist == "negbin"){
      model_formula <- as.formula(paste("read ~ treatment + offset(log(norm_factors))", if (!is.null(covariates)) paste("+", colnames(covariates), collapse = " + "), sep = " "))
      model = glm.nb(formula = model_formula, data = genewise_dataframe)
    }

    coefficients <- coef(summary(model))[, "Estimate"]
    p_values <- coef(summary(model))[, "Pr(>|z|)"]
    padj_values <- p.adjust(p_values, method = "BH")
    log2fold <- log(exp(abs(coefficients)), base = 2)

    c(coefficients, p_values, padj_values, log2fold)
  }

  stopImplicitCluster()

  colnames(result_GLM) <- c("Coefficient", "P_value", "P_adj_value", "Log2fold")

  if (padj){
    pvals <- result_GLM[, "P_adj_value"]
  } else {
    pvals <- result_GLM[, "P_value"]
  }

  log2fold <- result_GLM[, "Log2fold"]

  idx <- 1 * (pvals < pval_thre)

  if (l2fc){
    idx <- idx * (log2fold > l2fc_thre)
  }

  output <- list(DE_idx = idx, pvals = pvals, log2fc = log2fold)

  return(output)
}
