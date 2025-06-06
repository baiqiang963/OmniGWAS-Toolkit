% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDL.rg.parallel.R
\name{HDL.rg.parallel}
\alias{HDL.rg.parallel}
\title{High-definition likelihood inference of genetic correlations with multiple cores}
\usage{
HDL.rg.parallel(gwas1.df, gwas2.df, LD.path, Nref = 335265,
  N0 = min(gwas1.df$N, gwas2.df$N), output.file = "", numCores = 2,
  eigen.cut = "automatic", jackknife.df = FALSE, 
  intercept.output = FALSE, fill.missing.N = NULL, lim = exp(-18))
}
\arguments{
\item{gwas1.df}{A data frame including GWAS summary statistics of genetic variants for trait 1. 
The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS. If the GWAS is based on logistic regression, `b` should be the logarithm of OR (odds ratio) and `se` is the standard error of log(OR). Notice: SNPs with missing N will be removed.}

\item{gwas2.df}{A data frame including GWAS summary statistics of genetic variants for trait 2. 
The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS. If the GWAS is based on logistic regression, `b` should be the logarithm of OR (odds ratio) and `se` is the standard error of log(OR). Notice: SNPs with missing N will be removed.}

\item{LD.path}{Path to the directory where linkage disequilibrium (LD) information is stored.}

\item{Nref}{Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 335265.}

\item{N0}{Number of individuals included in both cohorts. The estimated genetic correlation is usually robust against misspecified N0. 
If not given, the default value is set to the minimum sample size across all SNPs in cohort 1 and cohort 2.}

\item{output.file}{Where the log and results should be written. If you do not specify a file, the log will be printed on the console.}

\item{numCores}{The number of cores to be used.}

\item{eigen.cut}{Which eigenvalues and eigenvectors in each LD score matrix should be used for HDL. 
Users are allowed to specify a numeric value between 0 and 1 for eigen.cut. For example, eigen.cut = 0.99 means using the leading eigenvalues explaining 99% of the variance
and their correspondent eigenvectors. If the default 'automatic' is used, the eigen.cut gives the most stable heritability estimates will be used.}

\item{jackknife.df}{Logical, FALSE by default. Should the block-jackknife estimates be returned?}

\item{intercept.output}{Logical, FALSE by default. Should the intercept terms be included in estimates.df?}

\item{fill.missing.N}{If NULL (default), the SNPs with missing N are removed. One of "median", "min" or "max" can be given so that the missing N will be filled accordingly. For example, "median" means the missing N are filled with the median N of the SNPs with available N.}

\item{lim}{Tolerance limitation, default is exp(-18)}

}

\value{
A list is returned with:
\itemize{
\item{rg }{The estimated genetic correlation.}
\item{rg.se }{The standard error of the estimated genetic correlation.}
\item{P }{P-value based on Wald test.}
\item{estimates.df }{A detailed matrix includes the estimates and standard errors of heritabilities, genetic covariance and genetic correlation.}
\item{eigen.use }{The eigen.cut used in computation.}
\item{jackknife.df }{Only if argument `jackknife.df` is true. A matrix includes the block-jackknife estimates of heritabilities, genetic covariance and genetic correlation.}
}
}
\description{
This function is a parallel version of HDL.rg. The function returns the estimate and standard error of the genetic correlation between two traits based on GWAS summary statistics. 
Note: The parallelization is implemented by doSNOW pakckage. Please make sure it has been successfully installed.
}
\note{
Users can download the precomputed eigenvalues and eigenvectors of LD correlation matrices for European ancestry population. The download link can be found at https://github.com/zhenin/HDL/wiki/Reference-panels
These are the LD matrices and their eigen-decomposition from 335,265 genomic British UK Biobank individuals. Three sets of reference panel are provided: 
1) 1,029,876 QCed UK Biobank imputed HapMap3 SNPs. The size is about 33 GB after unzipping. Although it takes more time, using the imputed panel provides more accurate estimates of genetic correlations. 
Therefore if the GWAS includes most of the HapMap3 SNPs, then we recommend using the imputed reference panel.
2) 769,306 QCed UK Biobank imputed HapMap2 SNPs. The size is about 18 GB after unzipping.If one of your GWAS includes most of the HapMap 2 SNPs, but many SNPs (more than 1%) in the above HapMap 3 reference panel are absent, 
then this HapMap2 panel is more proper to be used for HDL. 
3) 307,519 QCed UK Biobank Axiom Array SNPs. The size is about 7.5 GB after unzipping.
}
\examples{
\dontrun{
## The GWAS summary statistics for birth weight 
data(gwas1.example)

## The GWAS summary statistics for type 2 diabetes
data(gwas2.example)

## The path to the directory where linkage disequilibrium (LD) information is stored.
LD.path <- "/Users/zhengning/Work/HDL/package/UKB_array_SVD_eigen90_extraction"

res.HDL <- HDL.rg.parallel(gwas1.example, gwas2.example, LD.path, numCores = 2)
res.HDL
}
}
\references{
Ning Z, Pawitan Y, Shen X (2020). High-definition likelihood inference of genetic correlations 
across human complex traits. \emph{Accepted by Nature Genetics, waiting for publishing}.
}
\seealso{
HDL tutorial: https://github.com/zhenin/HDL
}
\author{
Zheng Ning, Xia Shen
}
