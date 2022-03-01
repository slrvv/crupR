#devtools::load_all("/project/wig/persia/CRUP_package/CRUP/R/")
context("Enhancer Prediction")
library(GenomicRanges)
library(rtracklayer)
##################################################################
# create input data
##################################################################

files.short <- c(system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))

inputs <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3)


#create the metaData frame
metaData.short <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
                       condition = c(2,2,2), replicate = c(1,1,1),
                       bamFile = files.short, inputFile = inputs)


pred.expected <- readRDS(file = system.file("extdata", "condition2_predictions.rds", package="crupR"))
#the actual normalized counts
data_matrix <- readRDS(file = system.file("extdata", "condition2_normalized.rds", package="crupR"))

#create a list
norm = list(metaData = metaData.short, D = data_matrix)

##################################################################
# test enhancerPrediction()
##################################################################

testthat::test_that("the error messages of getEnhancers() work", {

  testthat::expect_error(crupR::getEnhancers(data = norm, classifier = "/wrong/directory/classifier/", C = 4),
                         "/wrong/directory/classifier/ is not a valid directory")
})
testthat::test_that("getEnhancers() runs as expected",{
  pred <- crupR::getEnhancers(data = norm, C = 2)
  pred_short <- pred$D[which(as.character(GenomicRanges::seqnames(pred$D)) == "chr8")]
  testthat::expect_equal(length(pred), 2)
  testthat::expect_equal(pred_short$prob, pred.expected$prob, tolerance = 1e-5)
})

