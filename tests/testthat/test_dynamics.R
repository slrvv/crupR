#devtools::load_all("/project/wig/persia/CRUP_package/CRUP/R/")
context("Find condition-specific enhancers")
library(GenomicRanges)
library(rtracklayer)

##################################################################
# create input data
##################################################################

files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))

inputs <- rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3)

inputs2 <- rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3)

metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = c(inputs, inputs2))

metaData1 <- subset(metaData, condition == 1)
metaData2 <- subset(metaData, condition == 2)


pred1 <- readRDS(system.file("extdata", "condition1_predictions.rds", package = "crupR"))
#colnames(GenomicRanges::mcols(pred1)) <- c("prob")
pred1 <- list("metaData" = metaData1, "D" = pred1)
pred2 <- readRDS(system.file("extdata", "condition2_predictions.rds", package = "crupR"))
#colnames(GenomicRanges::mcols(pred2)) <- c("prob")
pred2 <- list("metaData" = metaData2, "D" = pred2)

preds <- list(pred1, pred2)

##################################################################
# test enhancerDynamics()
##################################################################

testthat::test_that("the error messages of getDynamics() work",{

  testthat::expect_error(crupR::getDynamics(data = preds, w_0 = 1.2, C = 10),
                         "1.2 is not in range [0,1].", fixed = TRUE)

  testthat::expect_error(crupR::getDynamics(data = preds, cutoff = -0.5, C = 10),
                         "-0.5 is not in range [0, 1].", fixed = TRUE)
})

testthat::test_that("getDynamics runs as expected",{
  dynamics.expected <- readRDS(system.file("extdata", "differential_enhancers.rds", package = "crupR"))
  dynamics <- crupR::getDynamics(data = preds, C = 2)
  testthat::expect_equal(length(dynamics), 2)
  testthat::expect_equal(length(dynamics$sumFile), 1)
  testthat::expect_equal(dynamics$sumFile$cond1_1, dynamics.expected$cond1_1, tolerance = 1e-9)
})



