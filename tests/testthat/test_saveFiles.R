#devtools::load_all("/project/wig/persia/CRUP_package/CRUP/R/")
context("Save the output files of crupR")
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


metaData2 <- subset(metaData, condition == 2)

#recreate normalize() output
data_matrix <- readRDS(file = system.file("extdata", "condition2_normalized.rds", package="crupR"))
norm = list("metaData" = metaData2, "D" = data_matrix)

#recreate getSE() output
pred2 <- readRDS(system.file("extdata", "condition2_predictions.rds", package = "crupR"))
peaks <- readRDS(file = system.file("extdata", "condition2_peaks.rds", package="crupR"))
cluster <- readRDS(file = system.file("extdata", "condition2_clusters.rds", package="crupR"))
pred2 <- list("metaData" = metaData2, "D" = pred2, "peaks" = peaks, "clusters" = cluster)

#recreate getDynamics output
sumFile <- readRDS(system.file("extdata", "differential_enhancers.rds", package = "crupR"))
dynamics <- list("metaData" = metaData, "sumFile" = sumFile)

#recreate getTargets output
units <- readRDS(system.file("extdata", "RegulatoryUnits.rds", package="crupR"))
targets <- list("metaData" = metaData, "Units" = units)

##################################################################
# test enhancerDynamics()
##################################################################

testthat::test_that("the error messages of saveFiles() work",{
  test.path <- paste0(system.file("extdata", package = "crupR"), "/")
  testthat::expect_error(crupR::saveFiles(data = pred2, modes = c("rds"), outdir = "/wrong/path/"),
                         "Directory /wrong/path/ doesn't exist!", fixed = TRUE)
  
  testthat::expect_error(crupR::saveFiles(data = pred2, modes = c("bigWigs"), outdir = test.path),
                         "One of the modes you chose is not supported by this function.\n Please choose from rds, bigWig, bed, bedGraph, beds and UCSC.", 
                         fixed = TRUE)
  
  testthat::expect_error(crupR::saveFiles(data = pred2, modes = c("UCSC"), outdir = test.path),
                         "Only outputs of getTargets() can be saved in the mode UCSC.", 
                         fixed = TRUE)
  
  testthat::expect_error(crupR::saveFiles(data = pred2, modes = c("beds"), outdir = test.path),
                         "Only outputs of getDynamics() can be saved in the mode beds.",
                         fixed = TRUE)
  
  testthat::expect_error(crupR::saveFiles(data = dynamics, modes = c("bedGraph"), outdir = test.path),
                         "Only outputs of getSE() can be saved in the modes bedGraph and/or bed.",
                         fixed = TRUE)
  
  testthat::expect_error(crupR::saveFiles(data = norm, modes = c("rds"), outdir = test.path),
                         "Outputs of normalize() can't be saved using this function, use saveRDS() instead.",
                         fixed = TRUE)
  
})

testthat::test_that("saveFiles runs as expected",{
  test.path <- paste0(system.file("extdata", package = "crupR"), "/")
  out.rds <- paste0(test.path, "prediction.rds")
  out.bw <- paste0(test.path, "prediction.bw")
  out.bedGraph <- paste0(test.path, "singleEnh.bedGraph")
  out.bed <- paste0(test.path, "clusterEnh.bed")
  out.beds <- paste0(test.path, "dynamicEnh__cluster_c1.bed")
  out.ucsc <- paste0(test.path, "RegulatoryUnits.interaction")
  files <- c(out.rds, out.bw, out.bedGraph, out.bed, out.beds, out.ucsc)
  
  crupR::saveFiles(data = pred2, modes = c("rds", "bigWig", "bedGraph", "bed"), outdir = test.path)
  crupR::saveFiles(data = dynamics, modes = c("beds"), outdir = test.path)
  crupR::saveFiles(data = targets, modes = c("UCSC"), outdir = test.path)
  
  test.file <- readRDS(out.rds)
  
  testthat::expect_true(all(file.exists(files)))
  testthat::expect_equal(pred2$D$prob, test.file$prob, tolerance = 1e-5)
  file.remove(files)
})



