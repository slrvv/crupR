#' Find the target genes of the regulatory, condition-specific clusters
#'
#' @param data Outputfile of the previous step (getDynamics())
#' @param expr GRanges object containing the gene expression counts of RNAseq experiments for each condition and its replicates
#' @param genome Genome used in the .bam files of the RNAseq experiments. Possible options are 'mm9', 'mm10', 'hg19' and 'hg38'.
#' @param TAD.file Path to the TAD file to use for finding the target genes. If set to NULL, the default file is used (only if the 'mm10' genome was used)
#' @param cutoff cutoff for correlation between cluster and gene. Default is 0.9.
#' @param nearest Logical: if set, the nearest gene is taken to build the regulatory regions.
#' @param C Number of cores to use
#' @return A list containing the meta data of the experiments
#' and a GRanges object containing the dynamic regulatory units
#' @examples
#' 
#' #first get the output of crupR::getDynamics so skip this
#' files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
#'            system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
#'            system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
#' inputs <- rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3)
#' inputs2 <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3)  
#' metaData <- data.frame(HM = rep(c("H3K4me1", "H3K4me3", "H3K27ac"),2),
#'                        condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
#'                        bamFile = files, inputFile = c(inputs, inputs2))
#' clusters <- readRDS(system.file("extdata", "differential_enhancers.rds", package="crupR"))
#' dynamics <- list("metaData" = metaData, "sumFile" = clusters)
#' #load your GRanges object containing the gene expressions counts and run the function
#' expr <- readRDS(system.file("extdata", "expressions.rds", package="crupR"))
#' getTargets(data = dynamics, expr = expr, genome = "mm10", C = 2)
#'
#' @export
#' @importFrom GenomicRanges mcols GRanges makeGRangesFromDataFrame start end promoters strand seqnames nearest distance
#' @importFrom GenomicFeatures genes
#' @import parallel
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom IRanges IRanges subsetByOverlaps %within%
#' @import BSgenome
#' @import TxDb.Mmusculus.UCSC.mm9.knownGene
#' @import TxDb.Mmusculus.UCSC.mm10.knownGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom utils read.table

getTargets <- function(data, expr, genome, TAD.file = NULL, cutoff = 0.9, nearest = FALSE, C = 1){
  start_time <- Sys.time()
  cat('\n')
  
  metaData <- data$metaData
  conds <- unique(metaData$condition)
  gr<- data$sumFile
  GenomeInfoDb::seqlevels(gr) <- paste0("chr", gsub("chr|Chr","", GenomeInfoDb::seqlevels(gr)))
  
  IDs <- list()
  for(i in seq_along(conds)){
    sub = subset(metaData, condition == conds[i])
    if(is.numeric(conds[i])) IDs[[i]] <- paste0("cond", conds[i], "_", unique(sub$replicate))
    else IDs[[i]] <- paste0(conds[i], "_", unique(sub$replicate))
  }
  
  if (!(genome %in% genome_values)) stop(paste0("Genome " , genome, " currently not provided. Choose one of:", paste(genome_values, collapse=',')));
  if(cutoff < 0.5 | cutoff > 1) stop(paste0(cutoff, " is not in range [0.5,1]."));

  if (genome == "mm10") {txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  } else if (genome == "mm9") {txdb  <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
  } else if (genome == "hg19") {txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if (genome == "hg38") {txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene	}

  if (nearest == FALSE){
      if (is.null(TAD.file) & genome == "mm10") {TAD.file   <- system.file("extdata", "mESC_mapq30_KR_all_TADs.bed", package = "crupR")
      } else if (is.null(TAD.file) & genome != 'mm10')  {stop(paste0("You have to provide your own file with TAD domains (fitting to the genome of choice)."))}
  }
	
  TAD <- GenomicRanges::GRanges()
  if (nearest == FALSE){
      TAD <- read.table(TAD.file, col.names = GR_header_short)
      TAD <- GenomicRanges::makeGRangesFromDataFrame(TAD[which((TAD$end-TAD$start) > 0),])
      GenomeInfoDb::seqlevels(TAD) = paste0("chr", gsub("chr|Chr","", GenomeInfoDb::seqlevels(TAD)))
  }
  units <- get_units(gr, expr, TAD, unlist(IDs), C, cutoff, txdb, nearest)
 
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return( list("metaData" = data$metaData,  "Units" = units))
}
