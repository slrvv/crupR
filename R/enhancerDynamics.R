#' Clusters the active enhancer into condition-specific groups
#'
#' @param data List containing the metadata of the experiments and
#' the bin-wise enhancer prediction values as a GRanges object 
#' (calculated in the prior step)
#' @param w_0 The minimum difference between the normalized prediction means
#' that two enhancers need in order to be included in the clustering. 
#' Default is 0.5 ([0,1]).
#' @param cutoff cutoff for the p-values calculated during the clustering.
#'  Default is 0.05 ([0,1]).
#' @param W Number of bins +/- the current bin that should be included 
#' when calculating the p-values. Default is 10 ([5, 50]).
#' @param C Number of cores to use
#' @return A list containing the meta data of the experiments and
#' the condition-specific clusters in a GRanges object
#' @examples                                            
#' #recreate the outputs of getEnhancers()
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
#' metaData1 <- subset(metaData, condition == 1)
#' metaData2 <- subset(metaData, condition == 2)
#' pred1 <- readRDS(system.file("extdata", "condition1_predictions.rds", package = "crupR"))
#' pred1 <- list("metaData" = metaData1, "D" = pred1)
#' pred2 <- readRDS(system.file("extdata", "condition2_predictions.rds", package = "crupR"))
#' pred2 <- list("metaData" = metaData2, "D" = pred2)
#' #put the outputs in a list
#' predictions <- list(pred1, pred2)
#' #run the function
#' getDynamics(data = predictions, C = 2)
#'
#' @export
#' @importFrom GenomicRanges GRanges elementMetadata merge mcols granges start end reduce findOverlaps seqnames makeGRangesFromDataFrame resize
#' @importFrom dplyr mutate lag filter setdiff
#' @importFrom matrixStats rowVars colMedians colMaxs
#' @import parallel
#' @importFrom IRanges values
#' @importFrom magrittr extract %>%
#' @importFrom stats ks.test
#' @importFrom utils combn getFromNamespace
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomeInfoDb keepSeqlevels

getDynamics <- function(data, w_0 = 0.5, cutoff = 0.05, W = 10, C = 1){
  start_time <- Sys.time()
  cat('\n')
  
  if(cutoff < 0 | cutoff > 1) stop(paste0(cutoff, " is not in range [0, 1]."))
  if(w_0 < 0 | w_0 > 1) stop(paste0(w_0, " is not in range [0,1]."))
  if(W < 10 | W > 50) stop(paste0('Choose a ', W, " in range [10,50]."))
  #if(!is.list(data) || ! unique(do.call(c, lapply(data, function(x) names(x)))) %in% c('metaData', 'D')) stop('Input must be a list of lists with names: \'metaData\' and \'D\'.')
  if(!is.list(data) | ! all(unique(do.call(c, lapply(data, function(x) names(x)))) %in% c('metaData', 'D'))) stop('Input must be a list of lists with names: \'metaData\' and \'D\'.')
  
  ##################################################################
  # read and combine probabilities:
  ##################################################################

  cat("Merge enhancer probabilites for all samples and conditions ...\n")
  
  metaData <- do.call(rbind, lapply(data, function(x) as.data.frame(x$metaData)))
  probsList <- lapply(data, function(x) {   tmp <- x$D
  							 GenomicRanges::mcols(tmp) <- GenomicRanges::mcols(tmp[,grep('prob$', colnames(mcols(tmp)))])
  							 return(tmp)
  							 })
  conds <- unique(metaData$condition)

  IDs <- list()
  for(i in seq_along(conds)){
    sub = subset(metaData, condition == conds[i])
    if(is.numeric(conds[i])) IDs[[i]] <- paste0("cond", conds[i], "_", unique(sub$replicate))
    else IDs[[i]] <- paste0(conds[i], "_", unique(sub$replicate))
  }
  IDs[[length(IDs)+1]] <- 'null'
  
  probs <- GenomicRanges::GRanges()
  for(i in seq_along(probsList)){
    m <- probsList[[i]]
    colnames(GenomicRanges::elementMetadata(m)) <- unlist(IDs)[i]
    probs <- GenomicRanges::merge(probs, m, all = TRUE)
  }
  probs$null <- 0
  
  ##################################################################
  # get pairwise p-values and call differential peaks
  ##################################################################

  cat("Calculate pairwise p-values ...\n")
  pvalues <- get_pairwisePvalues(p = probs, I = IDs, w_0 = w_0, W = W, p.thres = cutoff, C = C)
  
  cat("Get condition-specific enhancer peaks ...\n")
  probs <- get_cluster(p = probs, pvalues = pvalues, I = IDs)
  
  GenomicRanges::elementMetadata(probs)[,'null'] <- NULL
  IDs <- IDs[-length(IDs)]
  
  #if (length(unique(probs$cluster)) == 0  && unique(probs$cluster) == 0) stop(paste0("No significant peaks found between any two conditions.\n"))

  if (length(unique(probs$cluster)) == 0){
    stop(paste0("No significant peaks found between any two conditions.\n"))}

  cat("Combine peaks by significance pattern ...\n")
  peaks <- get_ranges(p = probs, I = IDs, W = W, C = C)
  
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(list('metaData' = metaData, 'sumFile' = peaks))

}
