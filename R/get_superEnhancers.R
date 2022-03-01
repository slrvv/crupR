#' Uses the enhancer predictions to find clusters of enhancers ('superenhancer' - like).
#'
#' @param data List containing the meta data of the experiments and the enhancer prediction values in a GRanges object calculated by enhancerPrediction()
#' @param cutoff Cutoff for enhancer probabilities. Default is 0.5.
#' @param distance Maximum distance (bp) for clustering. Default is 12500
#' @param C Number of cores to use
#' @return A list containing the meta data of the experiments, the enhancer prediction values as a GRanges object (like the input data),
#' the enhancer peak calls as a GRanges object (can be exported as a bedGraph) and the cluster of peaks in a GRanges object (can be exported as bed file)
#' @examples
#'
#' # first recreate the output of crupR:getPredictions
#' files <- c(system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
#'            system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
#' inputs <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3)
#' #create the metaData frame
#' metaData <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
#'                        condition = c(2,2,2), replicate = c(1,1,1),
#'                        bamFile = files, inputFile = inputs)
#' pred <- readRDS(system.file("extdata", "condition2_predictions.rds", package="crupR"))                                                                                           
#' prediction <- list(metaData = metaData, D = pred)
#' #run the function
#' getSE(data = prediction, C = 2)
#' @export
#' @importFrom GenomicRanges start reduce mcols findOverlaps width seqnames GRanges
#' @import parallel
#' @importFrom S4Vectors queryHits subjectHits

getSE <- function(data, cutoff = 0.5, distance = 12500, C = 1){
  start_time <- Sys.time()
  cat('\n')
  
  if(cutoff < 0 | cutoff > 1) stop(paste0(cutoff, " is not a valid cutoff. Please choose a cutoff between 0 and 1."))
  if(distance < 0) stop(paste0(distance, " is not a valid distance. Please choose a distance greater than 0."))
  #if(!is.list(data) || ! names(data) %in% c('metaData', 'D')) stop('Input must be a list with names: \'metaData\' and \'D\'.')
  if(!is.list(data) | ! all(names(data) %in% c('metaData', 'D'))) stop('Input must be a list with names: \'metaData\' and \'D\'.')
  
  peaks <- GenomicRanges::GRanges()
  peaks <- c(peaks, get_enhancerPeaks(data$D, cutoff, C))
  cat(paste0(length(peaks), " single enhancer peak(s) found.\n"))
  
  cluster <- GenomicRanges::GRanges()
  cluster <- c(cluster, get_enhancerCluster(peaks, distance, C))
  cat(paste0(length(cluster), " enhancer cluster(s) found.\n"))
  
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(list("metaData" = data$metaData, "D" = data$D, "peaks" = peaks, "cluster" = cluster))
}
