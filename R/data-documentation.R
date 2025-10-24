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
#' For "allcause" mortality data.frame contains age-specific person-year data, "male_py" and "female_py", in addition to "age", "male" and "female".
#'
#' @examples
#'  names(Mortality)       # WHO global regions
#'  names(Mortality[[1]])  # Sites for which baseline mortality rates are available
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
#' @details The list object for each risk model contains the 9 site-specific cancer mortality risk models derived from Life Span Study ("allsolid", "esophagus", "stomach", "colon", "liver", lung", "bladder", "breast", "leukaemia").
#'
#' Each site-specific data.frame contains information for the risk model (a vector of parameter estimates "para", a matrix of variance-covariance matrix of parameter estimates "var" and a function to calculate the risk "f").
#'
#' @examples
#'  names(LSS_mortality)   # Sites for which LSS mortality risk models are available
#'  names(LSS_mortality$allsolid)    # Available dose response models
#'  LSS_mortality$allsolid$L$err     # Linear ERR model for all solid cancer motality
#'  LSS_mortality$allsolid$L$ear     # Linear EAR model for all solid cancer motality
#'
#'  LSS_mortality$leukaemia$LQ$err   # Linear-quadratic ERR model for leukaemia motality
#'  LSS_mortality$lung$L$err         # Linear EAR model for lung cancer motality
#'
#'  # Plotting LSS all solid cancer mortality risk model
#'  plot_riskmodel( rm=LSS_mortality$allsolid$L, title="LSS all solid cancer mortality, Linear",  leg_pos=c(0.4, 0.95) )
#'  # Plotting LSS Leukaemia incidence risk model
#'  plot_riskmodel( rm=LSS_incidence$leukaemia$LQ, title="LSS leukaemia incidence", ymax=c(1.5, .3), add=c(0.01,0) )

#'
#' @seealso [LSS_incidence], [plot_riskmodel()]
#' @references Ozasa, K., Y. Shimizu, A. Suyama et al. Studies of the mortality of atomic bomb survivors, Report 14, 1950-2003: an overview of cancer and noncancer diseases. Radiat Res 177(3): 229-243 (2012).
#' @source IARC
"LSS_mortality"



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
#'
#'  names(Incidence)       # WHO global regions
#'  names(Incidence[[1]])  # Sites for which baseline incidence rates are available
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


#' LSS incidence risk models
#'
#' LSS incidence risk models for 13 cancer sites.
#'
#' @docType data
#' @usage LSS_incidence
#' @format A list object of :
#' \describe{
#'   \item{allsolid}{a list object which contains risk model information (see Details)}
#'   \item{leukaemia}{a list object which contains risk model information (see Details)}
#'   \item{esophagus}{a list object which contains risk model information (see Details)}
#'   \item{stomach}{a list object which contains risk model information (see Details)}
#'   \item{colon}{a list object which contains risk model information (see Details)}
#'   \item{liver}{a list object which contains risk model information (see Details)}
#'   \item{lung}{a list object which contains risk model information (see Details)}
#'   \item{prostate}{a list object which contains risk model information (see Details)}
#'   \item{pancreas}{a list object which contains risk model information (see Details)}
#'   \item{bladder}{a list object which contains risk model information (see Details)}
#'   \item{breast}{a list object which contains risk model information (see Details)}
#'   \item{thyroid}{a list object which contains risk model information (see Details)}
#'   \item{brainCNS}{a list object which contains risk model information (see Details)}
#' }
#' @details The list object for each risk model contains the 13 site-specific cancer incidence risk models derived from Life Span Study ("leukaemia", "esophagus", "stomach", "colon", "liver", "lung", "prostate", "pancreas", "bladder", "breast", "thyroid", "brainCNS").
#'
#' Each site-specific a list object contains information for the risk model (a vector of parameter estimates "para", a matrix of variance-covariance matrix of parameter estimates "var" and a function to calculate the risk "f").
#'
#' @examples
#'  names(LSS_incidence)   # Sites for which LSS incidence risk models are available
#'  names(LSS_incidence$allsolid)    # Available dose response models
#'  LSS_incidence$allsolid$L$err     # Linear ERR model for all solid cancer incidence
#'  LSS_incidence$allsolid$L$ear     # Linear EAR model for all solid cancer incidence
#'
#'  LSS_incidence$leukaemia$LQ$err   # Linear-quadratic ERR model for leukaemia incidence
#'  LSS_incidence$lung$L$err         # Linear ERR model for thyroid cancer incidence
#'
#'  # Plotting LSS all solid cancer mortality risk model
#'  plot_riskmodel( rm=LSS_incidence$allsolid$L, title="LSS all solid cancer mortality, Linear",  leg_pos=c(0.4, 0.95) )
#'
#'  # Plotting LSS Leukaemia incidence risk model
#'  plot_riskmodel( rm=LSS_incidence$leukaemia$LQ, title="LSS leukaemia incidence", ymax=c(1.5, .3), add=c(0.01,0) )
#'
#'  # Plotting LSS thyroid cancer incidence risk model
#'  plot_riskmodel( rm=LSS_incidence$thyroid$L, title="LSS thyroid cancer incidence", ymax=c(.5, .3) )

#'
#' @seealso [LSS_mortality], [plot_riskmodel()]
#' @references Ozasa, K., Y. Shimizu, A. Suyama et al. Studies of the mortality of atomic bomb survivors, Report 14, 1950-2003: an overview of cancer and noncancer diseases. Radiat Res 177(3): 229-243 (2012).
#' @source IARC
"LSS_incidence"









