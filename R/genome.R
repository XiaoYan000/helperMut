#  $$$$$$\  $$$$$$$$\ $$\   $$\  $$$$$$\  $$\      $$\ $$$$$$$$\
# $$  __$$\ $$  _____|$$$\  $$ |$$  __$$\ $$$\    $$$ |$$  _____|
# $$ /  \__|$$ |      $$$$\ $$ |$$ /  $$ |$$$$\  $$$$ |$$ |
# $$ |$$$$\ $$$$$\    $$ $$\$$ |$$ |  $$ |$$\$$\$$ $$ |$$$$$\
# $$ |\_$$ |$$  __|   $$ \$$$$ |$$ |  $$ |$$ \$$$  $$ |$$  __|
# $$ |  $$ |$$ |      $$ |\$$$ |$$ |  $$ |$$ |\$  /$$ |$$ |
# \$$$$$$  |$$$$$$$$\ $$ | \$$ | $$$$$$  |$$ | \_/ $$ |$$$$$$$$\
#  \______/ \________|\__|  \__| \______/ \__|     \__|\________|



#' Get chunks from regions
#'
#' From a genome file or a GR object generate a chunked version
#'
#' @param gr a GRanges object or a BSgenome object
#' @param wl the window length of the resulting chunks
#' @param unlist if result should be returned as a single GR, default true
#'
#' @return A GRanges object chunked
#' @export
#'
#' @examples
get_region_chunks <- function(gr,wl = 10000,unlist = TRUE) {

  # deciding which type of object we have
  if (is(gr,"GRanges")){
    start_vector = BiocGenerics::start(gr)
    end_vector = BiocGenerics::end(gr)
  } else if (is(gr,"BSgenome")) {
    start_vector = rep(1,length(gr))
    end_vector = GenomeInfoDb::seqlengths(gr)
  } else {
    stop("gr class not supported")
  }

  start = purrr::map2(.x = start_vector,
                      .y = end_vector,
                      .f = seq,
                      by = wl)

  end = purrr::map2(.x = start,
                    .y = end_vector,
                    function(x,max_value){
                      end = dplyr::lead(x) - 1
                      end[length(end)] = max_value
                      return(end)
                    })

  lol = list(
    as.character(GenomeInfoDb::seqnames(gr)),
    end,
    start
  )

  purrr::pmap(.l = lol,.f = function(name,end,start){
    GenomicRanges::GRanges(
      seqnames = as.character(name),
      ranges = IRanges::IRanges(start = start,end = end)
    )
  }) -> res

  if (unlist){
    names(res) = NULL
    # the warning s due to multiple things don't share seqnmas
    # could be solved by explicitely seting seqnames
    # I think it wouldn't be worth the effort.
    res = suppressWarnings(do.call(what = "c",args = res))
  }

  return(res)
}



