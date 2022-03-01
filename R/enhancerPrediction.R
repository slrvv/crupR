#' Predicts the occurence of enhancers using the normalized HMs counts.
#'
#' @param data List containing the metaData of the experiments that were normalized and also the normalized counts as GRanges object (basically the output of the normalize() step)
#' @param classifier The path of the classifier to use for the prediction. When set to NULL (default), the default clasifier is used.
#' @param all Output of separated probabilities. Default is FALSE.
#' @param C Number of cores to use. Default is 1.
#' @return A list containing the meta data of the experiments whose enhancers were predicted and
#' the enhancer probabilities for each 100 bp bin in the genome as a GRanges object
#' @examples
#' #first recreate the output of crupR::normalize (so skip this)
#' files <- c(system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
#'           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
#'           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
#' inputs <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR")) 
#' metaData <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
#'                     condition = c(2,2,2), replicate = c(1,1,1),
#'                     bamFile = files, inputFile = inputs)
#' d <- readRDS(system.file("extdata", "condition2_normalized.rds", package="crupR"))
#' norm <- list(metaData = metaData, D = d)
#' #let's run the actual function
#' getEnhancers(data = norm, C = 2)
#'
#' @export
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @import randomForest
#' @importFrom stats predict
#' @importFrom GenomicRanges mcols elementMetadata end reduce
#' @import parallel

getEnhancers <- function(data, classifier = NULL, all = FALSE,
                         C = 1){
  start_time <- Sys.time()
  cat('\n')
  
  if (!is.null(classifier) && (!dir.exists(classifier))) stop(paste0(classifier, " is not a valid directory"))
  #if(!is.list(data) || ! names(data) %in% c('metaData', 'D')) stop('Input must be a list with names: \'metaData\' and \'D\'.')
  if(!is.list(data) | !all(names(data) %in% c('metaData', 'D'))) stop('Input must be a list with names: \'metaData\' and \'D\'.')
  ##################################################################
  # Read classifier files and ecdf
  ##################################################################

  cat("Get classifier and empirical distribution function ...\n")

  if(is.null(classifier)){
  	classifierF1 <- system.file("extdata", "active_vs_inactive.rds", package = "crupR")
  	classifierF2 <- system.file("extdata", "enhancer_vs_active_promoter.rds", package = "crupR")
  } else{
  	classifierF1 <- file.path(classifier, "active_vs_inactive.rds")
  	classifierF2 <- file.path(classifier, "enhancer_vs_active_promoter.rds")
  }

  check_file(classifierF1)
  check_file(classifierF2)
  classifier1 <- readRDS(classifierF1)
  classifier2 <- readRDS(classifierF2)

  feat1 <- unique(gsub('_.*','',names(classifier1$forest$xlevels)))
  feat2 <- unique(gsub('_.*','',names(classifier2$forest$xlevels)))
  featAll <- unique(c(feat1, feat2))

  ecdf_file <- system.file("extdata", "ecdf.rds", package = "crupR")
  check_file(ecdf_file)
  ecdf <- readRDS(ecdf_file)

  ##################################################################
  # Quantile normalization
  ##################################################################

  cat(paste0("Quantile normalize counts for features ...\n"))
  
  d <-data$D
  dNorm <- data$D
  for (f in featAll) GenomicRanges::mcols(dNorm)[,f] <- preprocessCore::normalize.quantiles.use.target(matrix(GenomicRanges::mcols(dNorm)[,f]), get_targetQuantileNorm(ecdf[[f]]))

  ##################################################################
  # Extend data matrix
  ##################################################################

  cat("Extend data matrix ...\n")
  
  d_ext  <- extend_dataMatrix(N = 5, df = data.frame(d), f = feat1)
  zero.idx <- which(rowSums(d_ext[,-c(seq_len(3))]) == 0)
  dNorm_ext <- extend_dataMatrix(N = 5, df = data.frame(dNorm), f = featAll)
  dNorm_ext <- dNorm_ext[-zero.idx,]

  ##################################################################
  # make predictions
  ##################################################################

  cat("Predict enhancers for each 100 bp bin ...\n")
  
  mid <- round(nrow(dNorm_ext)/2, 0)
  pred1 <- predict(classifier1, dNorm_ext[seq_len(mid),], type = "prob")[,2]
  pred1 <- c(pred1, predict(classifier1, dNorm_ext[(mid+1):nrow(dNorm_ext),], type = "prob")[,2])
  pred2 <- predict(classifier2, dNorm_ext[seq_len(mid),], type = "prob")[,2]
  pred2 <- c(pred2, predict(classifier2, dNorm_ext[(mid+1):nrow(dNorm_ext),], type = "prob")[,2])

  GenomicRanges::elementMetadata(d) <- NULL
  GenomicRanges::mcols(d)["prob"] <- 0
  GenomicRanges::mcols(d[-zero.idx])["prob"] <- pred1 * pred2
  if(all==TRUE){
        GenomicRanges::mcols(d)["probA"] <- 0
        GenomicRanges::mcols(d)["probE"] <- 0
  	GenomicRanges::mcols(d[-zero.idx])["probA"] <- pred1
  	GenomicRanges::mcols(d[-zero.idx])["probE"] <- pred2
  }
  cat(paste0('time: ', format(Sys.time() - start_time), "\n"))
  return(list("metaData" = data$metaData ,"D" = d))
}
