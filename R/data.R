# $$$$$$$\   $$$$$$\ $$$$$$$$\  $$$$$$\
# $$  __$$\ $$  __$$\\__$$  __|$$  __$$\
# $$ |  $$ |$$ /  $$ |  $$ |   $$ /  $$ |
# $$ |  $$ |$$$$$$$$ |  $$ |   $$$$$$$$ |
# $$ |  $$ |$$  __$$ |  $$ |   $$  __$$ |
# $$ |  $$ |$$ |  $$ |  $$ |   $$ |  $$ |
# $$$$$$$  |$$ |  $$ |  $$ |   $$ |  $$ |
# \_______/ \__|  \__|  \__|   \__|  \__|


#' Ordered mutation types based on Cosmic signatures website
#'
#' A vector containing the ordered mutation types at k=1
#' from the sanger cosmic website
#'
#' @format A vector with 96 mutation types
#' @source \url{https://cancer.sanger.ac.uk/signatures_v2/Signature-2.png}
"order_ms96_cosmicSignatures"

#' A small VR object to use in examples.
#'
#' A small VRanges object to use in the examples of the documentation
#'
#' @format VRanges object
"vr_example"

#' A small GR object to use in examples.
#'
#' A small GRanges object to use in the examples of the documentation
#'
#' @format GRanges object
"gr_example"

#' IUPAC DNA codes
#'
#' A list of each DNA letter and its representative set in ACTG letters.
#'
#' @format A list with all IUPAC DNA codes and its meanings
#' @source \url{https://www.bioinformatics.org/sms/iupac.html}
"dna_codes"

#' Color scheme for CA mutations
#'
#' A dataset containing color codes for CA mutations.
#'
#' @format A vector or list of color codes
#' @usage data(tr_colors_CA)
#' @export
"tr_colors_CA"

#' Color scheme for CT mutations
#'
#' A dataset containing color codes for CT mutations.
#'
#' @format A vector or list of color codes
#' @usage data(tr_colors_CT)
#' @export
"tr_colors_CT"
