context("Find target genes")
library(GenomicRanges)
#library(data.table)

##################################################################
# create input data
##################################################################
#setwd("/project/wig/persia/CRUP_package/CRUP/tests/testthat/")
files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))

inputs <- c(rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3), rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3))

metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = inputs)

sumFile <- readRDS(system.file("extdata", "differential_enhancers.rds", package = "crupR"))
expr.gr <- readRDS(system.file("extdata", "expressions.rds", package = "crupR"))

dynamics <- list("metaData" = metaData, "sumFile" = sumFile)

##################################################################
# test enhancerTargets()
##################################################################


units_expected <- readRDS(system.file("extdata", "RegulatoryUnits.rds", package="crupR"))

testthat::test_that("the error messages of getTargets() work",{

  testthat::expect_error(crupR::getTargets(data = dynamics, expr = expr.gr, genome = "mm11", C = 6),
                         "Genome mm11 currently not provided. Choose one of:hg19,mm10,mm9,hg38")

  testthat::expect_error(crupR::getTargets(data = dynamics, expr = expr.gr, genome = "mm10", cutoff = 0.2, C = 6),
                         "0.2 is not in range [0.5,1].", fixed = TRUE)

  testthat::expect_error(crupR::getTargets(data = dynamics, expr = expr.gr, genome = "mm9", C = 6),
                         "You have to provide your own file with TAD domains (fitting to the genome of choice).", fixed = TRUE)
})

testthat::test_that("enhancerTargets runs as expected",{
  targets <- crupR::getTargets(data = dynamics, expr = expr.gr, genome = "mm10", C = 2)

  testthat::expect_equal(length(targets), 2)
  testthat::expect_equal(length(targets$Units), 2)
  testthat::expect_identical(colnames(GenomicRanges::mcols(targets$Units)), colnames(GenomicRanges::mcols(units_expected)))
  testthat::expect_identical(levels(GenomicRanges::mcols(targets$Units)$CORRELATED_GENE),
  			levels(GenomicRanges::mcols(units_expected)$CORRELATED_GENE))
})

ng_units_expected <- readRDS(system.file("extdata", "RegulatoryUnitsNearestGene.rds", package = "crupR"))

testthat::test_that("Nearest Gene mode runs as expected", {
  ng_targets <- crupR::getTargets(data = dynamics, expr = expr.gr, genome = "mm10", nearest = TRUE, C = 2)

  testthat::expect_equal(length(ng_targets), 2)
  testthat::expect_equal(length(ng_targets$Units), 1)
  testthat::expect_identical(colnames(mcols(ng_targets$Units)), colnames(mcols(ng_units_expected)))
  testthat::expect_identical(levels(GenomicRanges::mcols(ng_targets$Units)$NEAREST_GENE),
  			     levels(GenomicRanges::mcols(ng_units_expected)$NEAREST_GENE))
})

