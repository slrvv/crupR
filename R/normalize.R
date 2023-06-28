#' Input normalization of ChIPseq counts.
#'
#' @param metaData A data frame containing all important information about the ChIPseq experiments, i,e. for which histone modification they were conducted, path to the file, which condition or replicate, ...
#' @param condition The number of the condition to normalize.
#' @param replicate The number of the replicate to normalize.
#' @param genome Reference genome. Either a character (one of: mm9, mm10, hg19, hg38) or a Seqinfo object
#' @param mapq A number, Minimum mapping quality of the reads that should be included. Default is 10.
#' @param sequencing A string, what type of sequencing was used. Possible options: paired and single
#' @param input.free A boolean, Whether or not to run the function in the inputfree mode, in case there is no input experiment available. Default is False.
#' @param chroms A vector of strings. Specify the relevant chromosomes if needed. Then only the reads on these chromosomes are considered for normalization and used in further analysis. Default is 'all'.
#' @param C Number of cores to use. Default is 1.
#' @return A list containing the meta data of the experiments that were normalized and the normalized counts in a GRanges object
#' @examples
#'
#'files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
#'           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
#'           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
#'           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
#'           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
#'           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
#'
#'inputs <- c(rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3), 
#'            rep(system.file("extdata", "Condition2.Input.bam", package="crupR"),3))
#'                                                        
#'metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
#'                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
#'                       bamFile = files, inputFile = inputs)
#'                                             
#' normalize(metaData = metaData, condition = 1, replicate = 1,
#'     genome = "mm10", sequencing = "paired", C = 2)
#' 
#' #Example for a customized genome:
#' genome = GenomeInfoDb::Seqinfo(seqnames=c("chr3", "chr4", "chrM"),
#'               	seqlengths=c(1000, 2000, 500))
#'
#' @export
#' @importFrom bamsignals bamProfile
#' @importFrom Rsamtools scanBamHeader BamFile
#' @importFrom stats knots
#' @importFrom GenomicRanges seqinfo mcols tileGenome width
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevels
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
#' @importFrom BSgenome.Mmusculus.UCSC.mm9 BSgenome.Mmusculus.UCSC.mm9
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38

normalize <- function(metaData, condition, replicate, genome, mapq = 10, sequencing, input.free = FALSE, chroms = NULL, C = 1) {
  start_time <- Sys.time()
  cat('\n')
  #in general: comment more!!!!
  m <- metaData[which(metaData$condition == condition & metaData$replicate == replicate),]
  if(nrow(m) != 3){
    if(nrow(m) == 0) stop(paste0("The chosen combination of condition and replicate is not valid.\n There are no files for condition ", condition," replicate ", replicate))
    else if(nrow(m) > 3) stop(paste0("There are too many bam files for condition ", condition," replicate ", replicate))
    else stop(paste0("There are not enough bam files for condition ", condition," replicate ", replicate))
  }

  for(path in m$bamFile) check_file(path)
  if(! input.free) for(path in m$inputFile) check_file(path)
  
  hm <- c()
  for(i in seq_along(hm_values)) hm <- c(hm, normalizePath(as.vector(m[which(m$HM == hm_values[i]),]$bamFile)[1]))
  if(input.free == FALSE) for(i in seq_along(hm_values)) hm <- c(hm, normalizePath(as.vector(m[which(m$HM == hm_values[i]),]$inputFile)[1]))

  bamHM <- hm[seq_len(3)]
  if(input.free == TRUE){bamInp <- NULL
  }  else{bamInp <- unique(hm[4:length(hm)])}
  
  if (!is(genome, 'Seqinfo') && genome == "mm10") {genome <- GenomeInfoDb::seqinfo(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  }  else if (!is(genome, 'Seqinfo') && genome == "mm9") {genome  <- GenomeInfoDb::seqinfo(BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9)
  }  else if (!is(genome, 'Seqinfo') && genome == "hg19") {genome <- GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  }  else if (!is(genome, 'Seqinfo') && genome == "hg38") {genome <- GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  }  else if(!is(genome, 'Seqinfo')) {stop(paste0("Your genome is neither one of ",paste0(genome_values, collapse=',')," nor is it a valid Seqinfo object."))}

  ##################################################################
  # get and adjust binned genome
  ##################################################################

  cat("Prepare the binned genome ...\n")
  
  gr <- get_binned_genome(genome, chr=chroms)
  bf <- Rsamtools::BamFile(bamHM[1])
  si <- GenomicRanges::seqinfo(bf)
  if (!("chr1" %in% names(Rsamtools::scanBamHeader(bamHM[1])[[1]]$targets))){ 
    GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(si)[1]
  }
  
  ##################################################################
  # get counts from ChIP-seq experiments
  ##################################################################

  cat("Get summarized counts from ChIP-seq experiments ...\n")

  names <- c(hm_values, paste0("Input_",hm_values))
  if(input.free) names <- hm_values
  if (length(unique(bamInp)) == 1) {
    bamInp <- bamInp[1]
    names = c(hm_values, "Input_All")
  }
  
  if (sequencing == "paired") {
      counts <- parallel::mcmapply(
                bamsignals::bamProfile, 
                bampath = c(bamHM, bamInp),
                MoreArgs = list(gr = gr, 
                binsize = 100,
                mapqual = mapq,
                ss = FALSE,
                paired.end = "midpoint",
                filteredFlag = 1024,
                verbose = FALSE),
                mc.cores = C, 
                SIMPLIFY = FALSE)

  } else if (sequencing == "single") {
      counts <- parallel::mcmapply(
              bamsignals::bamProfile, 
              bampath = c(bamHM, bamInp),
              MoreArgs = list(gr = gr, 
              binsize = 100,
              mapqual = mapq,
              ss = FALSE,
              shift = 100,
              filteredFlag = 1024,
              verbose = FALSE),
              mc.cores = C, 
              SIMPLIFY = FALSE)
  } else {
    stop("Sequencing parameter is not valid.\n Choose one of:",
         paste(sequencing_values, collapse=','))
  }
  for(i in seq_along(counts)) counts[[i]] <- unlist(as.list(counts[[i]]))
  names(counts) <- names

  ##################################################################
  # Input normalizatrion of the ChIP-seq counts
  ##################################################################
  if(! input.free){
      cat("Normalize histone modifications by Input ...\n")

      if ("Input_All" %in% names(counts)) {
        countsNorm <- lapply(hm_values, 
                             function(x) { # Put it in a function girl
                               log2((counts[[x]] + 1)/(counts[[paste0("Input_All")]] + 1})))
      } else { 
        countsNorm <- lapply(hm_values, 
                             function(x) log2((counts[[x]] + 1)/(counts[[paste0("Input_",x)]] + 1)))
  } else { 
    countsNorm <- counts
  }
  ##################################################################
  # create data matrix for all normalized ChIP-seq experiments
  ##################################################################
  cat("Create summarized data matrix ...\n")

  GenomicRanges::mcols(gr) <- matrix(unlist(countsNorm), 
                                      ncol = 3,
                                      byrow = FALSE,
                                      dimnames = list(NULL, hm_values))
                             
  GenomeInfoDb::seqlevels(gr) <- paste0('chr', gsub('chr|Chr','',GenomeInfoDb::seqlevels(gr)))

  n <- GenomicRanges::mcols(gr)[,"H3K4me1"] + abs(min(GenomicRanges::mcols(gr)[,"H3K4me1"])) + 1
  d <- GenomicRanges::mcols(gr)[,"H3K4me3"] + abs(min(GenomicRanges::mcols(gr)[,"H3K4me3"])) + 1
  GenomicRanges::mcols(gr)[,"ratio"] <- log2(n/d)
  
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(list(metaData = m, D = gr))
}
