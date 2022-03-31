
##################################################################
# definition: valid parameter values:
##################################################################

hm_values <- c('H3K4me1', 'H3K4me3', 'H3K27ac')
sequencing_values <- c('paired', 'single')
genome_values <- c('hg19', 'mm10', 'mm9', 'hg38')

##################################################################
#definition: standard GRanges/DataFrame header:
##################################################################

GR_header <- c("seqnames", "start","end","width")
GR_header_short <- c("seqnames", "start","end")
DF_header <- c("chr", "start","end")

##################################################################
# function: re-header
##################################################################

reheader_DF <- function(DF, header){
  colnames(DF)[seq_along(header)] <- header
  return(DF)
}

##################################################################
# function: check if file exists
##################################################################

check_file <- function(f){
  if (!(file.exists(f))) {
    stop(paste0("File ",f," does not exist.\n"));
  }
}

##################################################################
# function: check if outdir exists
##################################################################

check_outdir <- function(d, alt_d){
  if (is.null(d)) d <- dirname(alt_d)
  if (!dir.exists(d)) {
    cat(paste0("Output directory '",d,"'is not a valid directory. \n
              Output directory is set to ",dirname(alt_d)));
    d <-paste0(dirname(alt_d),"/")
  }
  return(d)
}


##################################################################
# function: partition genome into bins
##################################################################

get_binned_genome <- function(genome, chr=NULL){
  if(! is.null(chr)) {GenomeInfoDb::seqlevels(genome) <- chr
  } else {GenomeInfoDb::seqlevels(genome) <- GenomeInfoDb::seqlevels(genome)[grep("^chr[0-9]{,2}$|chrX$", GenomeInfoDb::seqlevels(genome))]}
  binned <- GenomicRanges::tileGenome(genome, tilewidth = 100, cut.last.tile.in.chrom = TRUE)
  return(binned[-which(GenomicRanges::width(binned) != 100)])
}


##################################################################
# function: quantile normalize with target
##################################################################

get_targetQuantileNorm <- function(ecdf){
  x.ref <- knots(ecdf)
  y.ref <- ecdf(x.ref)
  temp <- c(0, round(y.ref*26337756, digits = 0))
  return(unlist(lapply(seq_along(x.ref), function(x) rep(x.ref[x], temp[x + 1] - temp[x]))))
}



##################################################################
# function: create extended data matrix (plus/minus) bins
##################################################################

extend_dataMatrix <- function(N, df, f){
  N_regions <- nrow(df)
  f_ext <- NULL
  for (i in seq_len(N)) {
    f_ext <- c(f_ext, paste0(f, "_left", i ))
    f_ext <- c(f_ext, paste0(f, "_right", i ))
  }
  df_ext <- cbind(df[,c(GR_header_short, f)], matrix(0, nrow = N_regions, ncol = length(f_ext)))
  colnames(df_ext) <- c(DF_header, f, f_ext)
  region_ext <- (N + 1):(N_regions - N)
  for (i in seq_len(N)) {
    df_ext[region_ext, f_ext[((2*i - 2)*length(f) + 1):((2*i - 1)*length(f))]] <- df[region_ext - i, f]
    df_ext[region_ext, f_ext[((2*i - 1)*length(f) + 1):(2*i*length(f))]] <-  df[region_ext + i, f]
  }
  return(df_ext)
}


#############################################################################
# function: sort peak candidates and remove overlapping peaks
# used in peak_calling function
#############################################################################

sort_peaks <- function(peaks){
  peaks <- peaks[sort(GenomicRanges::mcols(peaks)$prob, decreasing = TRUE, index.return = TRUE)$ix]
  count <- 0
  while (length(peaks) > (count + 1)) {
    count <- count + 1
    overlap.to <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(query = peaks[count], subject = peaks))
    if (length(overlap.to) == 1) next
    delete.index <- sort(overlap.to, decreasing = FALSE)[-1]
    peaks <- peaks[-delete.index]
  }
  return(peaks)
}


#############################################################################
# function: call peaks (from probabilities)
#############################################################################

get_enhancerPeaks <- function(gr, cutoff, C){
  if(length(gr) > 0){
	  candidates <- GenomicRanges::reduce(gr[which(GenomicRanges::mcols(gr)$prob > cutoff)])
	  peaks <- gr[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr, candidates))]
	  GenomicRanges::start(peaks) <- GenomicRanges::start(peaks) - 500
	  GenomicRanges::width(peaks) <- 1100
	  out <- parallel::mclapply(split(peaks, GenomicRanges::seqnames(peaks)), sort_peaks, mc.cores = C)
	  return(do.call("c", unname(out)))
  }
}


##################################################################
# function: call superenhancer-like candidates from called peaks
##################################################################

get_enhancerCluster <- function(peaks, peak.gap, C){
  if(length(peaks) > 0){
	  peaksRed <- GenomicRanges::reduce(peaks, min.gapwidth = peak.gap)
	  cluster <- peaksRed[which(GenomicRanges::width(peaksRed) > 1100)]
	  if (length(cluster) > 0) {
	    sort.max <- unlist(parallel::mclapply(seq_along(cluster),
		                        function(x) max(peaks$prob[S4Vectors::subjectHits(GenomicRanges::findOverlaps(cluster[x], peaks))]),
		                        mc.cores = C))
	    cluster <- cluster[sort(sort.max, index.return = TRUE, decreasing = TRUE)$ix]
	  }
	  return(cluster)
  }
}


##################################################################
# function: compute Kolmogorow-Smirnow test statistic
##################################################################

get_KS.STATISTIC <- function(d, N){
  Nplus <- rep(1/N, N*2)
  Nminus <- -1/N
  d.sort <- sort.int(d, index.return = TRUE, method = "quick")
  z <- cumsum(replace(Nplus, d.sort$ix <= N, Nminus))
  return(max(abs(z[c(diff(d.sort$x) != 0, TRUE)])))
}


##################################################################
# function: get indices of regions
#whose mean difference is higher than w_0
##################################################################

get_idx <- function(c1m, c2m, w_0, W){
  z <- rle( abs(c1m - c2m) < w_0) %>%
    unclass() %>%
    as.data.frame() %>%
    dplyr::mutate(end = cumsum(lengths), start = c(1, dplyr::lag(end)[-1] + 1)) %>%
    magrittr::extract(c(1,2,4,3)) %>%
    dplyr::filter(values == TRUE) %>%
    dplyr::filter(lengths >= (W*2+1))
  return(dplyr::setdiff(seq(length(c1m)), unlist(apply(z[,c('start', 'end')], 1, function(x) seq(x[1], (x[2]-W+1))))))
}

##################################################################
# function: get pairwise p values
##################################################################

get_pairwisePvalues <- function(p, I, w_0, W, p.thres, C){
  N <- W*2 + 1
  KS.FACTOR <- sqrt(N * 0.5)
  comb <- combn(seq(length(I)),2)
  ret <- list()

  for(i in seq(ncol(comb))){
    p.value <- rep(NA, length(p))
    i1 <- unlist(I[comb[1,i]])
    i2 <- unlist(I[comb[2,i]])

    if(length(i1) > 1){
      c1m <- rowMeans(as.matrix(GenomicRanges::mcols(p)[,i1]))
      c1s <- (1-sqrt(matrixStats::rowVars(as.matrix(GenomicRanges::mcols(p)[,i1]))) + 0.0000001)
    }else{
      c1m <- GenomicRanges::mcols(p)[,i1]
      c1s <- rep(1, length(c1m))
    }
    if(length(i2) > 1){
      c2m <- rowMeans(as.matrix(GenomicRanges::mcols(p)[,i2]))
      c2s <- (1-sqrt(matrixStats::rowVars(as.matrix(GenomicRanges::mcols(p)[,i2]))) + 0.0000001)
    }else{
      c2m <- GenomicRanges::mcols(p)[,i2]
      c2s <- rep(1, length(c2m))
    }
    z <- cbind(c1m/c1s, c2m/c2s)
    
    idx <- get_idx(c1m, c2m, w_0, W)
    if(length(idx) <= (W*2+1)) next
    idx <- idx[(W+1):(length(idx)- W)]
    
    D <- unlist(parallel::mclapply( idx,
                                     FUN = function(x)  get_KS.STATISTIC( z[(x - W):(x + W),], N),
                                     mc.cores = C,
                                     mc.allow.recursive = FALSE)
    )
    C_pKS2 <- utils::getFromNamespace("C_pKS2", "stats")
    p.value[idx] <- 1 - .Call(C_pKS2, p = KS.FACTOR * D, tol = 0.000001)
    p.adjust <- p.adjust(p.value, method = 'bonferroni')
    idx.significant <- which(p.adjust <= p.thres)
    if(length(idx.significant) == 0) stop("No regions with significant p-values were found.")

    this.result <- GenomicRanges::granges(p)[idx.significant]
    GenomicRanges::mcols(this.result) <- data.frame(p.value = p.value[idx.significant],
                                                        	p.adj = p.adjust[idx.significant],
                                                        	p.direction = unlist(lapply(idx.significant, function(x) (mean(z[(x - W):(x + W),1]) - mean(z[(x - W):(x + W),2]) <= 0))),
                                                        	idx = idx.significant,
                                                        	comparison = paste0(unique(gsub('_.*','',unlist(I[comb[,i]]))), collapse = ',')
        )
    ret[[ paste0(comb[,i], collapse=',')]] <- this.result[which(this.result$p.direction == TRUE)]
    ret[[ paste0(rev(comb[,i]), collapse=',')]] <- this.result[which(this.result$p.direction == FALSE)]
  }
  return(ret)
}

##################################################################
# function: get the pairwise pattern
##################################################################

get_cluster <- function(p, pvalues, I){
  comb <- combn(seq(seq_along(I)),2)
  comb.f <- cbind(comb, comb[c(2,1),])
  comb.f <- comb.f[,-which(comb.f[2,]==max(comb.f))]
  idx <- unique(unlist(lapply(pvalues, function(x) x$idx)))
  m <- matrix(0, nrow = length(idx), ncol = ncol(comb.f))
  colnames(m) <- apply(comb.f, 2, function(x) paste(x, collapse = ","))
  rownames(m) <- idx
  for (i in colnames(m)) m[as.character(pvalues[[i]]$idx),i] <- 1
  m=cbind(m[, which(comb.f[1,] != max(comb.f))], apply(m[, which(comb.f[1,] == max(comb.f))],1,prod))
  p$sign <- FALSE
  p$sign[idx] <- TRUE
  p$cluster <- 0
  p$cluster[idx] <- apply(m[,-max(comb.f)], 1, function(x) x%*%(2^(seq(length(x))-1))) 
  p$cluster[intersect(which(p$cluster == 0), idx[which(m[,max(comb.f)] == 1)])] <- 'ubiquitous'
  return(p)
}

##################################################################
# function: get summarized ranges
##################################################################

get_ranges <- function(p, I, W, C){
 	p.split <- split(p, GenomicRanges::seqnames(p), drop = TRUE)
	rle <- lapply(p.split, function(x) rle(x$cluster))
	d <- data.frame( seqnames=rep(names(rle), unlist(lapply(rle, function(x) length(x$lengths)))),
				end.idx=unlist(lapply(rle, function(x) cumsum(x$lengths))))
	d$start.idx <- d$end.idx - unlist(lapply(rle, function(x) x$lengths)) + 1
	d$cluster <- unlist(lapply(rle, function(x) x$values))
        
    #for(i in seq_along(p.split)){
    #    c <- levels(d$seqnames)
    #    p.split[[i]] <- GenomeInfoDb::keepSeqlevels(p.split[[i]],c, pruning.mode="tidy")}
    
	d$start <- unlist(lapply(p.split, function(x) GenomicRanges::start(x)[d$start.idx[which(d$seqnames ==unique(GenomicRanges::seqnames(x)))]]))
	d$end <- unlist(lapply(p.split, function(x) GenomicRanges::end(x)[d$end.idx[which(d$seqnames ==unique(GenomicRanges::seqnames(x)))]]))
	gr <- GenomicRanges::makeGRangesFromDataFrame(d[,-grep('.idx', colnames(d))], keep.extra.columns=TRUE)
	
	ovrlp <- GenomicRanges::findOverlaps(p, gr)[which(p$cluster !=0)]
	ovrlp.gr<-p[S4Vectors::queryHits(ovrlp)]
	ovrlp.split <- S4Vectors::split(ovrlp.gr, S4Vectors::subjectHits(ovrlp))
	names(ovrlp.split) <- gr$cluster[unique(S4Vectors::subjectHits(ovrlp))]
	ovrlp.red <- GenomicRanges::reduce(ovrlp.split)
	ovrlp.ext <- GenomicRanges::resize(ovrlp.red, width = width(ovrlp.red)+(W*100*2), fix = "center")
	ovrlp.final <- findOverlaps(p, ovrlp.ext)
	split.final <- S4Vectors::split(p[queryHits(ovrlp.final)], subjectHits(ovrlp.final))
	
	ret <- unlist(ovrlp.red)
	ret$cluster <- names(ovrlp.red)
	tab <- rev(table(ret$cluster))
	pattern.names <- paste0('c',seq(length(tab)))
	names(pattern.names) <- names(tab)
	ret$cluster <- pattern.names[ret$cluster]
	
	GenomicRanges::mcols(ret)[, unlist(I)] <- do.call("rbind",
                                                      parallel::mclapply(
                                                      split.final,
                                                      FUN = function(x) matrixStats::colMaxs(as.matrix(GenomicRanges::mcols(x)[,unlist(I)])),
                                                      mc.cores = C,
                                                      mc.allow.recursive = FALSE
                                                      )
  	)
	return(ret)
}

##################################################################
# function: generate colors
##################################################################
get_colors <- function(n,s=0.8,v=0.7){#,seed=11) {
  #set.seed(seed)
  h <- runif(1)
  H <- vector("numeric",n)
  for(i in seq_len(n)) {
    h <- (h + 0.618033988749895) %% 1
    H[i] <- h
  }
  hsv(H,s=s,v=v)
}

##################################################################
# function: plot summary of the K first differential enhancer regions
##################################################################
#' plots boxplots of the median enhancer prediction values of the enhancers in the condition-specific clusters
#'
#' @param D The getDynamics() output file of containing the GRanges object with the differential enhaners
#' @param num_plots Maximal number of cluster whose plots should be displayed (clusters are sorted by their sizes). This parameter can be set in case the number of clusters is very high (20+). Default is 20.
#' @return a ggplot2 object containing a boxplot for each cluster
#' @examples
#' #get the output of getDynamics()
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
#' clusters <- readRDS(system.file("extdata", "differential_enhancers.rds", package="crupR"))
#' dynamics <- list("metaData" = metaData, "sumFile" = clusters)
#' plotSummary(dynamics)
#'
#' @export
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot theme_bw guides xlab ylab
#' @importFrom reshape2 melt
#' @importFrom GenomicRanges mcols
#' @importFrom grDevices hsv
#' @importFrom stats runif

plotSummary <- function(D, num_plots = 9){
  p <- D$sumFile
  I <- unique(D$metaData$condition)
  R <- unique(D$metaData$replicate)
  colors <- get_colors(length(I))
    
  if(length(unique(p$cluster)) > num_plots)  p <- p[which(p$cluster %in% paste0("c",seq_len(num_plots)))]
  
  d <- suppressMessages(reshape2::melt(data.frame(GenomicRanges::mcols(p))))
  tab <- table(p$cluster)
  d$condition <- factor(paste0('cond ',gsub('_.*|cond','',d$variable)), paste0('cond ',levels=I))
  d$replicate <- factor(gsub('.*_','',d$variable), levels=R)
  d$label <- paste0(p$cluster, ' (', tab[p$cluster], ' regions)') 

  
  p <- ggplot2::ggplot(d, ggplot2::aes(x=condition, y=value, col=condition, shape=replicate)) +
    ggplot2::facet_wrap(label ~ . , ncol=round(length(unique(d$label))/4+0.3))+
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(	axis.title.x=ggplot2::element_blank(),
        			axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5),
        			axis.ticks.x=ggplot2::element_blank()) + 
    ggplot2::guides(shape='none', col='none') +
    ggplot2::xlab('') +
    ggplot2::ylab('max. probability') +
    ggplot2::scale_color_manual(values=colors, guide = ggplot2::guide_legend())
  
  return(p)
}


##################################################################
# function: create TAD if not defined
##################################################################

check_TAD <- function(t, TAD, this.region, regions){
  if (length(t) == 0) {
    precedeT <- GenomicRanges::precede(TAD, this.region)
    followT  <- GenomicRanges::follow(TAD, this.region)

    if(length(which(!is.na(precedeT))) == 0 ){start <- 1
    }else{	start <- GenomicRanges::end(TAD[rev(which(!is.na(precedeT)))[1]])
        	if (length(start) == 0) start <- 1}

    if(length(which(!is.na(followT))) == 0 ){
        end <- max(GenomicRanges::end(regions[which(GenomicRanges::seqnames(regions) == GenomicRanges::seqnames(this.region))]))
    }else{
        end <- GenomicRanges::start(TAD[(which(!is.na(followT)))[1]])
        if (length(end) == 0) end <- max(GenomicRanges::end(regions[which(GenomicRanges::seqnames(regions) == GenomicRanges::seqnames(this.region))]))}
    t <- GenomicRanges::GRanges(GenomicRanges::seqnames(this.region), IRanges::IRanges(start, width = (end - start + 1)))
  }
  return(t)
}

##################################################################
# function: correlate probabilities and gene exression counts
##################################################################

get_correlation <- function(i, cutoff, regions.gr, expr.gr, TAD.gr, I){
  interactions <- data.frame(stringsAsFactors = FALSE)
  this.region <- regions.gr[i]
  this.TAD <- IRanges::subsetByOverlaps(TAD.gr, this.region)
  this.TAD <- check_TAD(this.TAD, TAD.gr, this.region, regions.gr)
  this.genes.idx <- IRanges::'%within%'(expr.gr, this.TAD)
  
  if (sum(this.genes.idx) > 0) {
    this.genes <- expr.gr[this.genes.idx,]
    cor <- apply(as.matrix(GenomicRanges::mcols(this.genes)[,I]), 1,
                 function(x) cor(x, as.numeric(unlist(GenomicRanges::mcols(this.region)[,I]))))

    for (c in seq_along(cor)) {
      if (!is.na(cor[c]) && cor[c] >= cutoff) {

      promoter_start <- GenomicRanges::start(GenomicRanges::promoters(this.genes[c]))
      promoter_end <- GenomicRanges::end(GenomicRanges::promoters(this.genes[c]))
      if(as.character(GenomicRanges::strand(this.genes[c]))=="-"){
          promoter_start <- GenomicRanges::end(GenomicRanges::promoters(this.genes[c]))
          promoter_end <- GenomicRanges::start(GenomicRanges::promoters(this.genes[c])) }

      if(length(this.TAD)!= length(this.genes[c])){
        gene.TAD <- IRanges::subsetByOverlaps(TAD.gr, this.genes[c])
        gene.TAD <- check_TAD(gene.TAD, TAD.gr, this.genes[c], this.genes)
      } else { gene.TAD <- this.TAD }
      interactions <- rbind(interactions,
                            data.frame(data.frame(this.region)[,c(GR_header_short, "cluster", I)],
                                                  TAD_COORDINATES = paste0(gene.TAD),
                                                  CORRELATED_GENE = paste(GenomicRanges::mcols(this.genes)[c,"gene_id"]),
                                                  CORRELATED_GENE_CHR = GenomicRanges::seqnames(this.genes[c]),
                                                  CORRELATED_GENE_PROMOTER_START = promoter_start,
                                                  CORRELATED_GENE_PROMOTER_END = promoter_end,
                                                  CORRELATION = cor[c] ))
      }
    }
    if (length(interactions) > 0) return(GenomicRanges::makeGRangesFromDataFrame(interactions, keep.extra.columns = TRUE))
  }
}


##################################################################
# function: get nearest gene for a region
##################################################################

get_nearest_gene <- function(i, regions.gr, genes, I){
  this.region <- regions.gr[i]
  nearest <- genes[GenomicRanges::nearest(this.region, genes)]
  distance.to.nearest <- GenomicRanges::distance(this.region, nearest)
  promoter_start <- GenomicRanges::start(GenomicRanges::promoters(nearest))
  promoter_end <- GenomicRanges::end(GenomicRanges::promoters(nearest))
  if(as.character(GenomicRanges::strand(nearest))=="-"){
      promoter_start <- GenomicRanges::end(GenomicRanges::promoters(nearest))
      promoter_end <- GenomicRanges::start(GenomicRanges::promoters(nearest))
  }
  interactions <- data.frame(data.frame(this.region)[,c(GR_header_short, "cluster", I)],
                              NEAREST_GENE = GenomicRanges::mcols(nearest)[,"gene_id"],
                              NEAREST_GENE_CHR = GenomicRanges::seqnames(nearest),
                              NEAREST_GENE_PROMOTER_START = promoter_start,
                              NEAREST_GENE_PROMOTER_END = promoter_end,
                              DISTANCE_TO_NEAREST = distance.to.nearest
                              )
  return(GenomicRanges::makeGRangesFromDataFrame(interactions, keep.extra.columns = TRUE))
}


##################################################################
# function: correlate probabilities with gene expression values
#(in same TAD; per cluster)
##################################################################

get_units <- function(regions.gr, expr.gr, TAD.gr, I, C, cutoff, txdb, nearest = FALSE){
  if(nearest == FALSE){
      list <- parallel::mclapply(seq(length(regions.gr)),
      function(x) get_correlation(x, cutoff, regions.gr, expr.gr, TAD.gr, I),
      mc.cores = C)
  }else{
      suppressMessages(genes <- GenomicFeatures::genes(txdb))
      GenomeInfoDb::seqlevels(genes, pruning.mode = "coarse") = GenomeInfoDb::seqlevels(regions.gr)
      GenomeInfoDb::seqlengths(regions.gr) <- GenomeInfoDb::seqlengths(genes)
      list <- parallel::mclapply(seq(length(regions.gr)),
                                  function(x) get_nearest_gene(x, regions.gr, genes, I),
                                  mc.cores = C)
  }
  units <-  do.call("c", unname(unlist(list)))
  return(units)
}
