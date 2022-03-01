context("Normalization of the ChIPseq counts")
library(GenomicRanges)
##################################################################
# create input data
##################################################################

files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))

inputs <- c(rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3), rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3))

files.wrong <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           "/wrong/directory/ConditionX.bam",
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))

#create the metaData frame
metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = inputs)

metaData.wrong <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                             condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                             bamFile = files.wrong, inputFile = inputs) #one path is not right

metaData.if <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                          condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                          bamFile = files)#input experiment are missing

##################################################################
# test normalize()
##################################################################

#load the results to compare
norm.expected <- readRDS(file = system.file("extdata", "condition2_normalized.rds", package="crupR"))
norm.if.expected <- readRDS(file = system.file("extdata", "condition2_inputfree_normalized.rds", package="crupR"))

testthat::test_that("the error messages of normalize() work", {
  testthat::expect_error(crupR::normalize(metaData = metaData, condition = 2, replicate = 2, genome = "mm10", sequencing = "paired", C = 10),
                 "The chosen combination of condition and replicate is not valid.\n There are no files for condition 2 replicate 2")

  testthat::expect_error(crupR::normalize(metaData = metaData, condition = 2, replicate = 1, genome = "mm10", sequencing = "singles", C = 10),
                  "Sequencing parameter is not valid.\n Choose one of:paired,single")

  testthat::expect_error(crupR::normalize(metaData = metaData, condition = 2, replicate = 1, genome = "mm11", sequencing = "paired", C = 10),
                  "Your genome is neither one of hg19,mm10,mm9,hg38 nor is it a valid Seqinfo object.")

  testthat::expect_error(crupR::normalize(metaData = metaData.wrong, condition = 2, replicate = 1, genome = "mm10", sequencing = "paired", C = 10),
                         "File /wrong/directory/ConditionX.bam does not exist.")
})

testthat::test_that("The output file is okay",{
  norm <- crupR::normalize(metaData = metaData, condition = 2, replicate = 1, genome = "mm10", sequencing = "paired", chroms = c("chr8"), C = 2)
  norm_short <- norm$D[which(as.character(GenomicRanges::seqnames(norm$D))=="chr8")]
  testthat::expect_equal(length(norm), 2)
  testthat::expect_identical(norm$metaData, subset(metaData, condition == 2))
  testthat::expect_equal(norm_short$ratio, norm.expected$ratio)
  testthat::expect_equal(norm_short$H3K27ac, norm.expected$H3K27ac)
  testthat::expect_equal(norm_short$H3K4me3, norm.expected$H3K4me3)
  testthat::expect_equal(norm_short$H3K4me1, norm.expected$H3K4me1)
})

testthat::test_that("the inputfree mode runs w/o mistakes",{
  #run crupR::normalize in the non inputFree mode
  norm.if <- crupR::normalize(metaData = metaData, condition = 2, replicate = 1, genome = "mm10", sequencing = "paired", input.free = TRUE, chroms = c("chr8"), C = 2)
  norm.if_short <- norm.if$D[which(as.character(GenomicRanges::seqnames(norm.if$D))=="chr8")]
  testthat::expect_equal(length(norm.if), 2)
  testthat::expect_identical(norm.if$metaData, subset(metaData, condition == 2))
  testthat::expect_equal(norm.if_short$ratio, norm.if.expected$ratio)
})
