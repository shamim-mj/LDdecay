#' Example Genotype Dataset
#'
#' A small soybean genotype dataset included with LDdecay for
#' demonstrating linkage disequilibrium analysis and visualization.
#'
#' @format A data frame containing example soybean lines and SNP markers.
#' The \code{taxa} column identifies the genotype, and the remaining
#' columns contain SNP genotype dosages coded as 0, 1, or 2.
#'
#' @examples
#' data(geno)
#' head(geno)
"geno"


#' Example Physical Marker Map
#'
#' A small soybean physical marker map corresponding to the SNP markers
#' included in the \code{geno} example dataset.
#'
#' @format A data frame with three columns:
#' \describe{
#'   \item{SNP}{SNP marker identifier.}
#'   \item{Chromosome}{Chromosome identifier.}
#'   \item{Position}{Physical position in base pairs.}
#' }
#'
#' @examples
#' data(map)
#' head(map)
"map"