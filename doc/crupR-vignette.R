## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>")

## -----------------------------------------------------------------------------
files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
                     
inputs <- c(rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3), rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3))

## -----------------------------------------------------------------------------
metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = inputs)

metaData

## -----------------------------------------------------------------------------
normalized_1 <- crupR::normalize(metaData = metaData, condition = 1, replicate = 1, genome = "mm10", 
                                sequencing = "paired", C = 2)
normalized_2 <- crupR::normalize(metaData = metaData, condition = 2, replicate = 1, genome = "mm10", 
                                sequencing = "paired", C = 2)

## -----------------------------------------------------------------------------
normalized_1$metaData
normalized_1$D

## -----------------------------------------------------------------------------
normalized_1_inputFree <- crupR::normalize(metaData = metaData, condition = 1, replicate = 1,
                                          genome = "mm10", sequencing = "paired", input.free = TRUE, C = 2)
normalized_1_inputFree$metaData
normalized_1_inputFree$D

## -----------------------------------------------------------------------------
normalized_1_chr8 <- crupR::normalize(metaData = metaData, condition = 1, replicate = 1, genome = "mm10", 
                                sequencing = "paired", chroms = c("chr8"), C = 2)

## -----------------------------------------------------------------------------
prediction_1 <- crupR::getEnhancers(data = normalized_1, C = 2)
prediction_2 <- crupR::getEnhancers(data = normalized_2, C = 2)

## ---- eval = FALSE------------------------------------------------------------
#  prediction_1_own_class <- crupR::getEnhancers(data = normalized_1, classifier = "path/to/classifier", C = 2)

## -----------------------------------------------------------------------------
se <- crupR::getSE(data = prediction_2, C = 2)
se_strict <- crupR::getSE(data = prediction_2, cutoff = 0.7, C = 2)
se_close <- crupR::getSE(data = prediction_2, distance=10000, C = 2)

## -----------------------------------------------------------------------------
predictions <- list(prediction_1, prediction_2)

## -----------------------------------------------------------------------------
clusters <- crupR::getDynamics(data = predictions, C = 2)

## -----------------------------------------------------------------------------
#meta data
clusters$metaData

#clusters
clusters$sumFile

## -----------------------------------------------------------------------------
crupR::plotSummary(clusters)

## -----------------------------------------------------------------------------
expression <- readRDS(file = system.file("extdata", "expressions.rds", package="crupR"))

## -----------------------------------------------------------------------------
expression

## -----------------------------------------------------------------------------
targets <- crupR::getTargets(data=clusters, expr = expression, genome = "mm10", C = 2)

## -----------------------------------------------------------------------------
path_to_bed <- system.file("extdata", "mESC_mapq30_KR_all_TADs.bed", package="crupR")
targets <- crupR::getTargets(data = clusters, expr = expression, genome = "mm10", TAD.file = path_to_bed, C = 2)

## -----------------------------------------------------------------------------
targets_nearest <- crupR::getTargets(data=clusters, expr = expr, genome = "mm10", nearest = TRUE, C = 2)

## -----------------------------------------------------------------------------
targets$Units

## ---- eval = FALSE------------------------------------------------------------
#  out_dir <- ""
#  #save the GRanges object of the getEnhancers() step
#  saveFiles(data = prediction_1, modes = c("bigWig", "rds"), outdir = out_dir)
#  #save the GRanges object of the getSE() step
#  saveFiles(data = se, modes = c("bedGraph", "bed"), outdir = out_dir)
#  #save the GRanges object of the getDynamics() step
#  saveFiles(data = clusters, modes = "beds", outdir = out_dir)
#  #save the GRanges object of the getTargets() step
#  saveFiles(data = targets, modes = "UCSC", outdir = out_dir)

## ---- eval = FALSE------------------------------------------------------------
#  saveFiles(data = targets_nearest, modes = "UCSC", outdir = out_dir, nearest = TRUE)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

