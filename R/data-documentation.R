#' Cancer mortality rates of WHO global regions
#'
#' Cancer mortality rates of WHO global regions.
#'
#' @docType data
#' @usage Mortality
#' @format A list object of the five WHO global regions:
#' \describe{
#'   \item{"Aus-NZ Europe Northern America"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Northern Africa - Western Asia"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Latin America and Caribbean"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Asia excl. Western Asia"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Sub-Saharan Africa"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#' }
#' @details The list object for each region contains the 18 site-specific cancer mortality rates ("esophagus", "stomach", "colon", "liver", "pancreas", "lung", "breast", "prostate", "bladder", "brainCNS", "thyroid", "all_leukaemia", "all_cancer", "allsolid-NMSC", "allsolid", "leukaemia", "allcause", "survival").
#'
#' Each site-specific data.frame contains variables "age", "male" and "female".
#'
#' For "allcause" mortality data.frame contains person-year data "male_py" and "female_py" in addition to "age", "male" and "female".
#'
#' @examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
#'  # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence)
#'
#'  # Example 1: All solid cancer mortality rates for Region-1
#'  Mortality[[1]]$allsolid
#'
#'  # Example 2: Leukaemia mortality rates for Region-3
#'  Mortality[[3]]$leukaemia
#'
#'  # Example 3: A;;ll;-cause mortality rates for Region-5
#'  Mortality[[5]]$allcause
#'
#'  # Example 4: plotting lung cancer mortality rates
#'  plot_refdata( dat=Mortality, outcome="lung", leg_pos=c(0.27,0.95) )
#'
#' @seealso [Incidence], [plot_refdata()]
#' @source IARC
"Mortality"

#' Cancer incidence rates of WHO global regions
#'
#' Cancer incidence rates of WHO global regions.
#'
#' @docType data
#' @usage Incidence
#' @format A list object of the five WHO global regions:
#' \describe{
#'   \item{"Aus-NZ Europe Northern America"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Northern Africa - Western Asia"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Latin America and Caribbean"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Asia excl. Western Asia"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#'   \item{"Sub-Saharan Africa"}{a list object of 18 items, each of which contains a data.frame of cancer rates (see Details)}
#' }
#' @details The list object for each region contains the 18 site-specific cancer mortality rates ("esophagus", "stomach", "colon", "liver", "pancreas", "lung", "breast", "prostate", "bladder", "brainCNS", "thyroid", "all_leukaemia", "all_cancer", "allsolid-NMSC", "allsolid", "leukaemia").
#'
#' Each site-specific data.frame contains variables "age", "male" and "female".
#'
#'
#' @examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
#'  # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence)
#'
#'  # Example 1: All solid cancer incidence rates for Region-1
#'  Incidence[[1]]$allsolid
#'
#'  # Example 2: Leukaemia incidence rates for Region-3
#'  Incidence[[3]]$leukaemia
#'
#'  # Example 3: A;;ll;-cause incidence rates for Region-5
#'  Incidence[[5]]$allcause
#'
#'  # Example 4: plotting lung cancer incidence rates
#'  plot_refdata( dat=Incidence, outcome="lung", leg_pos=c(0.27,0.95) )
#'
#' @seealso [Mortality], [plot_refdata()]
#' @source IARC
"Incidence"


#' LSS mortality risk models
#'
#' LSS mortality risk models for 9 cancer types.
#'
#' @docType data
#' @usage LSS_mortality
#' @format A list object of :
#' \describe{
#'   \item{allsolid}{a list object which contains risk model information (see Details)}
#'   \item{esophagus}{a list object which contains risk model information (see Details)}
#'   \item{stomach}{a list object which contains risk model information (see Details)}
#'   \item{colon}{a list object which contains risk model information (see Details)}
#'   \item{liver}{a list object which contains risk model information (see Details)}
#'   \item{lung}{a list object which contains risk model information (see Details)}
#'   \item{bladder}{a list object which contains risk model information (see Details)}
#'   \item{breast}{a list object which contains risk model information (see Details)}
#'   \item{leukaemia}{a list object which contains risk model information (see Details)}
#' }
#' @details The list object for each region contains the 18 site-specific cancer mortality rates ("esophagus", "stomach", "colon", "liver", "pancreas", "lung", "breast", "prostate", "bladder", "brainCNS", "thyroid", "all_leukaemia", "all_cancer", "allsolid-NMSC", "allsolid", "leukaemia", "allcause", "survival").
#'
#' Each site-specific data.frame contains variables "age", "male" and "female".
#'
#' For "allcause" mortality data.frame contains person-year data "male_py" and "female_py" in addition to "age", "male" and "female".
#'
#' @examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
#'  # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence)
#'
#'  # Example 1: All solid cancer mortality rates for Region-1
#'  Mortality[[1]]$allsolid
#'
#'  # Example 2: Leukaemia mortality rates for Region-3
#'  Mortality[[3]]$leukaemia
#'
#'  # Example 3: A;;ll;-cause mortality rates for Region-5
#'  Mortality[[5]]$allcause
#'
#'  # Example 4: plotting lung cancer mortality rates
#'  plot_refdata( dat=Mortality, outcome="lung", leg_pos=c(0.27,0.95) )
#'
#' @seealso [Incidence], [plot_refdata()]
#' @references Ozasa, K., Y. Shimizu, A. Suyama et al. Studies of the mortality of atomic bomb survivors, Report 14, 1950-2003: an overview of cancer and noncancer diseases. Radiat Res 177(3): 229-243 (2012).
#' @source IARC
"LSS_mortality"








