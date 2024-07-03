# $$\   $$\ $$$$$$$$\ $$$$$$\ $$\       $$$$$$\
# $$ |  $$ |\__$$  __|\_$$  _|$$ |     $$  __$$\
# $$ |  $$ |   $$ |     $$ |  $$ |     $$ /  \__|
# $$ |  $$ |   $$ |     $$ |  $$ |     \$$$$$$\
# $$ |  $$ |   $$ |     $$ |  $$ |      \____$$\
# $$ |  $$ |   $$ |     $$ |  $$ |     $$\   $$ |
# \$$$$$$  |   $$ |   $$$$$$\ $$$$$$$$\\$$$$$$  |
#  \______/    \__|   \______|\________|\______/




.lenunique <- function(x){
  length(unique(x))
}

# STATS ===============


#' Get k-mer frequency in mutation-neighbouring regions
#'
#' From mutations in a VRanges object, obtain the k-mer oligonucleotide
#' frequencies at a given window length.
#'
#' @param vr a VRanges object with SNVs
#' @param k the extension of mutation context in one direction
#' @param wl the window length from which each mutation would be expanded
#' @param gr a GRanges object delimiting the region of interest
#' @param genome a BSgenome object
#' @param ... arguments for oligonucleotideFrequency
#'
#' @return
#' A integer named vector with the oligonucleotides selected at a given k.
#' @export
#'
#' @examples
#' get_k_freq(vr_example,
#'           k = 1,
#'           wl = 1000,
#'           genome = genome_selector("Hsapiens.1000genomes.hs37d5"))
get_k_freq <- function(vr,k,wl,genome,...){

  up = round(wl/2,0)
  down = round(wl/2,0)
  extended_region = extend(x = vr,upstream = up, downstream = down)

  # we reduce to avoid overlaps
  extended_region = GenomicRanges::reduce(extended_region)
  GenomeInfoDb::seqlevels(extended_region) = GenomeInfoDb::seqlevels(genome)
  gL = GenomeInfoDb::seqlengths(genome)
  GenomeInfoDb::seqlengths(extended_region) <- gL
  extended_region = GenomicRanges::trim(extended_region)

  # we now compute the frequency
  seqN = BSgenome::getSeq(genome, extended_region)
  K = (k*2)+1
  tri_freqs = Biostrings::oligonucleotideFrequency(seqN, width=K,...)
  tri_counts = apply(tri_freqs,2, sum)
  return(tri_counts)
}

#' @describeIn get_k_freq Get k-mer frequency in regions
get_k_freq_fromRegion <- function(k,gr,...){

  # this is not really needed but just in case
  extended_region = extend(x = gr,upstream = k, downstream = k)
  extended_region = GenomicRanges::reduce(extended_region)
  GenomeInfoDb::seqlevels(extended_region) = GenomeInfoDb::seqlevels(genome)
  gL = GenomeInfoDb::seqlengths(genome)
  GenomeInfoDb::seqlengths(extended_region) <- gL
  extended_region = GenomicRanges::trim(extended_region)

  # we now compute the frequency
  seqN = BSgenome::getSeq(genome, extended_region)
  K = (k*2)+1
  tri_freqs = Biostrings::oligonucleotideFrequency(seqN, width=K,...)
  tri_counts = apply(tri_freqs,2, sum)
  return(tri_counts)
}
