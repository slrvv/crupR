#' Function to save the different output objects of each step
#'
#' @param data An output list of one of the steps 2 to 4 or the get_superEnhancer step
#' @param modes Formats in which the GRanges object should be saved. Following modes are available:
#' for the output of getEnhancers(): "bigWig" for a bigWig file, "rds" for an .rds file
#' for the output of getDynamics(): "beds" for saving each cluster in a seperate bed file
#' for the output of getTargets(): "UCSC" for a UCSC interaction file
#' for the output of getSE(): "bedGraph" for saving the peak calls in a bedGraph file, "bed" for saving the clusters of peaks in a bed file
#' @param outdir Output directory in which the files should be saved
#' @param nearest Only relevant, if you want to save the output of enhancerTargets. Specifies if the output was produced by the nearest gene mode of the function or not. Default is false.
#' @return Nothing
#' @examples
#' \dontrun{
#' saveFiles(data = prediction_1_1, modes = c("bigWig", "rds"), outdir = "/example/dir/")
#' saveFiles(data = dynamics, modes = "beds", outdir ="/example/dir/")
#' saveFiles(data = targets, modes = "UCSC", outdir = "/example/dir/", nearest = TRUE)}
#' @export
#' @importFrom rtracklayer export
#' @importFrom GenomicRanges mcols seqnames
#' @importFrom utils write.table


saveFiles <- function(data, modes, outdir, nearest = FALSE){
  possible.modes <- c("rds", "bigWig", "bed", "bedGraph", "beds", "UCSC")
  modes.pred = c("bigWig","bedGraph", "rds","bed")
  if(!dir.exists(outdir)) stop(paste0("Directory ", outdir, " doesn't exist!"))
  if(!all(modes %in% possible.modes)){
      stop("One of the modes you chose is not supported by this function.\n Please choose from rds, bigWig, bed, bedGraph, beds and UCSC.")} 
  if(length(data) < 4 & any(modes %in% c("bedGraph", "bed"))) stop("Only outputs of getSE() can be saved in the modes bedGraph and/or bed.")
  if(!is.null(data$D)){
    if(all(colnames(GenomicRanges::mcols(data$D)) != "prob")){
      stop("Outputs of normalize() can't be saved using this function, use saveRDS() instead.")
    }
  } 
  if(is.null(data$sumFile) & any(modes %in% c("beds"))) stop("Only outputs of getDynamics() can be saved in the mode beds.")
  if(is.null(data$Units) & any(modes %in% c("UCSC"))) stop("Only outputs of getTargets() can be saved in the mode UCSC.")
  
  # for: enhancerPrediction
      if("bedGraph" %in% modes) {
          out.bedgraph <- paste0(outdir, "singleEnh.bedGraph")
          rtracklayer::export(data$peaks, out.bedgraph)
       }
       if("bed" %in% modes) {
           out.bed <- paste0(outdir, "clusterEnh.bed")
           rtracklayer::export(data$cluster, out.bed)
       }
       if("bigWig" %in% modes) {
	         data.D <- data$D
           if(length(GenomicRanges::mcols(data.D)) > 1){
             GenomicRanges::mcols(data.D) <- NULL
             metaCols <- GenomicRanges::mcols(data$D)
             #for the normal predictios
             out.bw <- paste0(outdir,"prediction.bw")
             data.D$score <- metaCols$prob
             rtracklayer::export(data.D, out.bw)
             #for probA
             out.bw <- paste0(outdir,"prediction_probA.bw")
             data.D$score <- metaCols$probA
             rtracklayer::export(data.D, out.bw)
             #for probE
             out.bw <- paste0(outdir,"prediction_probE.bw")
             data.D$score <- metaCols$probE
             rtracklayer::export(data.D, out.bw)
          } else {
            colnames(GenomicRanges::mcols(data.D)) <- c("score")
            out.bw <- paste0(outdir, "prediction.bw")
            rtracklayer::export(data.D, out.bw)}
       }
       if("rds" %in% modes){
           out.rds <- paste0(outdir, "prediction.rds")
           saveRDS(data$D, out.rds)
       }


  # for: enhancerDynamics
      if("beds" %in% modes){
        peaks <- data$sumFile
        clusters <- unique(peaks$cluster)

        for(c in clusters){
            if(! grepl("r", c)){
            out.bed <- paste0(outdir, paste0("dynamicEnh__cluster_",c,".bed"))
            write.table(data.frame(peaks)[which(peaks$cluster == c), GR_header_short], file = out.bed, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
            }
        }
      }

  #for: enhancerTargets
     if("UCSC" %in% modes){
       out.interaction <- ""
         units <- data$Units
         header <- "track type=interact name=\"Dynamic Promoter-Enhancer Pairs\" description=\"Dynamic Promoter-Enhancer Pairs\" interactDirectional=true visibility=full"
         enhancer <- as.matrix(data.frame(units)[,c("start","end")])
         promoter <-c()
         score <- 0
         if (nearest == FALSE){
             out.interaction <- paste0(outdir, paste0("RegulatoryUnits.interaction"))
             promoter <- as.matrix(data.frame(GenomicRanges::mcols(units)$CORRELATED_GENE_PROMOTER_START, GenomicRanges::mcols(units)$CORRELATED_GENE_PROMOTER_END))
             score <- GenomicRanges::mcols(units)$CORRELATION
         }else{
             out.interaction <- paste0(outdir, paste0("RegulatoryUnitsNearestGene.interaction"))
             promoter <- as.matrix(data.frame(GenomicRanges::mcols(units)$NEAREST_GENE_PROMOTER_START, GenomicRanges::mcols(units)$NEAREST_GENE_PROMOTER_END))
             score <- GenomicRanges::mcols(units)$DISTANCE_TO_NEAREST
         }

         interaction <- cbind(as.character(GenomicRanges::seqnames(units)),
                                                 apply(cbind(enhancer,promoter),1,min),
                                                 apply(cbind(enhancer,promoter),1,max),
                                                 rep(".", length(units)),
                                                 rep(0, length(units)),
                                                 score,
                                                 rep(".", length(units)),
                                                 rep("#7A67EE", length(units)),
                                                 as.character(GenomicRanges::seqnames(units)),
                                                 enhancer,
                                                 rep(".", length(units)),
                                                 rep(".", length(units)),
                                                 as.character(GenomicRanges::seqnames(units)),
                                                 promoter,
                                                 rep(".", length(units)),
                                                 rep(".", length(units)) )

        write.table(header, file = out.interaction, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
        write.table(interaction, file = out.interaction, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
     }

}
