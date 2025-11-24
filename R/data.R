#' GOTeDNA Example Data
#'
#' A subset of pre-formatted GOTeDNA metabarcoding data for testing and demonstration
#' of the GOTeDNA R package functions outside of the GOTeDNA Shiny application.
#'
#' @format ## 'D_mb_ex'
#' A data frame with 231 rows and 23 columns:
#' \describe{
#' \item{protocol_ID}{Common project ID to be used when lab and qPCR methods are identical, even for different species and locations}
#' \item{protocolVersion}{To be entered when project lab or field methods differ slightly}
#' \item{samp_name}{Sample Identifier from original field sample}
#' \item{eventID}{An identifier for the set of information associated with Sampling Event}
#' \item{primer}{Primer set used for amplification}
#' \item{species}{Target species name}
#' \item{domain}{Target species taxonomy}
#' \item{kingdom}{Target species kingdom}
#' \item{phylum}{Target species phylum}
#' \item{class}{Target species class}
#' \item{order}{Target species order}
#' \item{family}{Target species family}
#' \item{genus}{Target species genus}
#' \item{organismQuantity}{Number of DNA sequences reads in the sample}
#' \item{date}{Date of sample collection}
#' \item{LClabel}{When applicable (i.e., data collection on/by Indigenous peoples/lands), Local Contexts label for respecting and asserting Indigenous data sovereignty under the principles of OCAP (see: https://localcontexts.org/labels/biocultural-labels/). Under review.}
#' \item{decimalLongitude}{Longitude of sampling station}
#' \item{decimalLatitude}{Latitude of sampling station}
#' \item{station}{Sampling station identifier}
#' \item{year}{Year of sample collection}
#' \item{month}{Month of sample collection}
#' \item{detected}{Numeric, binary detected/not detected (1/0) based on minimum sequence copy threshold (msct)}
#' \item{msct}{Logical, where minimum sequence copy threshold = 10. FALSE filtered prior to analysis}
#' }
#' @source Fisheries and Oceans Canada; contact Anais.Lacoursiere@dfo-mpo.gc.ca for details.
"D_mb_ex"
