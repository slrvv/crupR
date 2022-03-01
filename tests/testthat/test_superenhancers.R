context("Find super enhancers")
library(GenomicRanges)
library(rtracklayer)

##################################################################
# create input data
##################################################################

files.short <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
                 system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
                 system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"))

inputs <- rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3)


#create the metaData frame
metaData.short <- data.frame(HM = c("H3K4me1","H3K4me3","H3K27ac"),
                             condition = c(1,1,1), replicate = c(1,1,1),
                             bamFile = files.short, inputFile = inputs)

pred <- readRDS(file = system.file("extdata", "condition2_predictions.rds", package="crupR"))
colnames(GenomicRanges::mcols(pred)) <- c("prob")
prediction <- list("metaData" = metaData.short, "D" = pred)

##################################################################
# test get_superEnhancers()
##################################################################

testthat::test_that("the error messages of getSE() work", {
  testthat::expect_error(crupR::getSE(data = prediction, cutoff = -1),
                         "-1 is not a valid cutoff. Please choose a cutoff between 0 and 1.")
  testthat::expect_error(crupR::getSE(data = prediction, distance = -100),
                        "-100 is not a valid distance. Please choose a distance greater than 0.")
})

#se <- crupR::getSE(data = prediction, C = 2)
testthat::test_that("getSE() runs as expected",{
  peaks.expected <- readRDS(file = system.file("extdata", "condition2_peaks.rds", package="crupR"))
  cluster.expected <- readRDS(file = system.file("extdata", "condition2_clusters.rds", package="crupR"))
  se <- crupR::getSE(data = prediction, C = 2)
  testthat::expect_equal(length(se), 4)
  testthat::expect_equal(length(se$peaks), 2)
  testthat::expect_equal(length(se$cluster), 1)
  testthat::expect_identical(GenomicRanges::start(se$peaks), GenomicRanges::start(peaks.expected))
  testthat::expect_identical(GenomicRanges::start(se$cluster), GenomicRanges::start(cluster.expected))
})
