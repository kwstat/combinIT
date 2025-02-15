#'
#' Copy number variation (CNV).
#'
#' This data set are about copy number variation (CNV) between normal and tumor tissue samples among six dogs.
#' In this data set, the value of CNV was measured as a signal intensity obtained from a comparative genomic hybridization (CGH) array,
#' with higher signals corresponding to higher copy number; see Franck et al. (2013) and Franck and Osborne (2016).
#' The data set was selected from 5899 sets (the full data have been made available as the supplementary material of the paper published by Franck et al. (2013)).
#' The test of interaction between the dogs and tisuues is of interest.
#'
#' @name CNV
#' @format A matrix with six rows (Dogs) and two columns (Tissues):
#' \describe{
#'   \item{Row1}{Dog1}
#'   \item{Row2}{Dog2}
#'   \item{Row3}{Dog3}
#'   \item{Row4}{Dog4}
#'   \item{Row5}{Dog5}
#'   \item{Row6}{Dog6}
#'   \item{Column1}{Normal tissue}
#'   \item{Column2}{Tumor}
#' }
#' @references
#' \enumerate{
#' \item
#'   Franck, C., Nielsen, D., Osborne, J.A. (2013). A method for detecting hidden additivity in two-factor unreplicated experiments.
#'    Computational Statistics and Data Analysis 67:95-104.
#' \item
#'  Franck, C., Osborne, J.A. (2016).  Exploring Interaction Effects in Two-Factor Studies using the hidden Package in R.
#'    R Journal 8 (1):159-172.
#' }
#'
NULL
