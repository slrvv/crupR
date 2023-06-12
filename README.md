# **crupR - Short Description**   

crupR is improved version of the enhancer detection pipeline CRUP ([C]ondition-specific [R]egulatory [U]nits [P]rediction). Its workflow consists of the same three main steps as the original pipeline (enhancer prediction, condition-specific enhancer detection, target gene detection) and an additional pre-preapring step (normalization). It uses different layers of epigenetic information in the form of ChIP-seq data from histone modification to provide a comprehensive list of regulatory units
consisting of dynamically changing enhancers and target genes.

#### *Installation*
Install crupR using following commad:

```
devtools::install_git("https://github.com/akbariomgba/crupR")
```


#### *Contact*

omgba@molgen.mpg.de

#### *Citation*

Since the crupR application note has not been submitted yet, you can cite this repository and the original paper if you use crupR:

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1860-7



## *crupR Overview*

*crupR* is a pipeline which aims to efficiently combine all its functions into a simple workflow. The main steps are:

* `normalize()`: Normalizes the ChIPseq counts
* `getPrediction()`: Uses the normalized ChIpseq counts to predict the occurence of active enhancers by applying a random forest based classifier
* `getDynamics()`: Uses the predicted probabilities to find dynamically changing, condition-specific enhancer clusters
* `getTargets(()`: Uses RNAseq experiments to find possible target genes for the dynamic enhancers.

Beyond that, *crupR* offers an additional function to summarize the single enhancers into super enhancers (`getSE()`), a function to visualise the condition-specific enhancer clusters (`plotSummary()`) and a function to export all the files produced by each *crupR* step as appropriate formats (i.e. bed, BigWig or bedGraph).

For a more detailed overview of *crupR* with more code examples, please check the vignette.

### *Getting started*

#### What will I need?

In order to run *crupR*, you need ChIPseq experiments for the histone modifications H3K27ac, H3K4me1 and H3K4me3. Additionally, it is recommended to also inclue an input experiment for those experiments, however if these aren't available for some reason, you can run *crupR* without them. All the experiments need to be already aligned and indexed.

#### Prepare the metaData file

After you prepared your ChIPseq files, you need to prepare the meta data frame. This data frame contains the paths to all the necessary bam files. It also contains further information about them, such as condition, replicate, histone modification and input file. 

Let's look at the metaData frame for the example files.
First, we'll need the paths to the bam files and input files.

```{r}
files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
                     
inputs <- c(rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3), rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3))
```

Now we can build the meta data frame by adding the information about the hinstone modification, condition and replicate for each file. As there exists only one sample per condition, all of the files get replicate "1".

```{r}
metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = inputs)

metaData
```
After creating the meta data frame, you are ready to run *crupR*.


### *Run the crupR pipeline*

#### *Step 0: Normalize the ChIPseq counts*

*crupR* needs normalized ChIPseq counts for its prediction, so we offer an additional function to normalize the ChIPseq files beforehand. You need to normalize each ChIPseq files for each replicate of each condition, which means we need to run *crupR* twice for our example files. Please note, that *crupR* can only normalize bam files for the genomes mm9, mm10, hg19 and hg38 or for a costum genome (which needs to be a Seinfo objects).  

```{r}
normalized_1 <- crupR::normalize(metaData = metaData, condition = 1, replicate = 1, genome = "mm10", 
                                sequencing = "paired", C = 2)
normalized_2 <- crupR::normalize(metaData = metaData, condition = 2, replicate = 1, genome = "mm10", 
                                sequencing = "paired", C = 2)
```
The output of `normalize()`  should be a list of length 2 containing the meta data of the samples that were normalized and a GRanges object containing the normalized ChIPseq counts for the binned genome of your choice. 
Per default, the files are input normalized. However, in case input experiments were not conducted or input files are not available for other reasons, *crupR* also offers the possibility of an input free mode. The output file will have the same structure.
Additionally, *crupR* offers the possibility to specify certain chromosmes that should be considered during the normalization, while discarding the remaining chromosomes. This might be relevant for users who wish to look for enhancers on specific chromosomes or who use BAM files that contain reads on only a few selected chromosomes. Specifying the chromsomes can improve runtime and also yields a smaller GRanges objects.
To specify the chromosomes the parameter *chroms* needs to be set by using a vector that contains the relevant chrosome names. You need to make sure, that the style of the chromosome names matches the style of the chromsome names in the respective bam file. I.e. if a BAḾ file uses the prefix "chr", then the chromosome names in the vector need to include the prefix, too. The resulting GRanges objects can be used as input for the next step. However, the subsequent steps will only consider the information on the specified chromosome(s).

#### *Step 1: Predict active enhancers with crupR*

Now, the normalized ChIPseq counts can be used to predict the occurrence of active enhancers on the genome of your choice.
Again, you'll have to run this step for each replicate of each condition, basically for each output of the normalization step.

```{r}
prediction_1 <- crupR::getEnhancers(data = normalized_1, C = 2)
prediction_2 <- crupR::getEnhancers(data = normalized_2, C = 2)
```

Per default, this function uses a classifier consisting of two random forests to predict the probability of an bin being an active enhancer. This default classifier was trained on mESC data. It is also possible to use your own classifier if you wish to do so. In that case you must specify the directory containing the classifiers. However, you need to make sure that the directory contains two classifiers that predict the same events that the default random forests predict. That means one classifer has to classify active vs. inactive genomic ranges, while the other one must classify enhancer vs active promoter regions. Additionally, the two classifiers also need to be named "active_vs_inactive.rds" and "enhancer_vs_active_promoter.rds" respectively. This is necessary so that *crupR* is able to accurately identify them.   

The output of this step is a list containing the truncated meta data file containing the information about the respective condition and replicate and a GRanges file containing the the binned genome with the enhancer prediction values for each bin. 

#### *Step 1.5: Find enhancer peaks and super enhancers with crupR*

*crupR* offers the `getSE()` function as an additional step after the enhancer prediction. Enhancers that are in close proximity are summarized into super enhancers during this step. However, this step does not solely return these super enhancers, but also all enhancer peak calls.

`getSE()` uses the output of the last step (enhancerPrediction), meaning the list with the meta data and the GRanges object, as input. There are two additnal input parameters to specifify the process:
* "cutoff": a threshold for the prediction values of the peaks (default: 0.5). The higher the cutoff, the stricter the peak calling process.
* "distance": the maximimum distance (bp) between peaks for clustering into superenhancers (default:12500). The lower the distance, the closer the single enhancers or peaks need to be in order to be clustered together.
A short example:
```{r}
se <- crupR::getSE(data = prediction_2, C = 2)
```
The output of this step is a list containing the same meta data file as the input, the same GRanges object that was produced by `getPrediction()`, a GRanges object containing the single enhancer peak calls (can be saved as a .bedGraph file) and a GRanges object containing the super enhancers or rather peak clusters (can be saved as a bed file).

This step is not necessary for the next ones and can be skipped if one isn't interested in the peak calls or peak clusters.

#### *Step 2: Find conditon-specific enhancer clusters*

In the second workflow step, the *crupR* function `getDynamics()` defines dynamic enhancer regions by applying a Kolmogorov-Smirnov test directly on the enhancer probabilities in a pairwise manner.

Before running the actual function, a small preparation step is necessary. The output lists of the previous step are used as input, but need to be contained by one file. Thus, they must all be put into one list.
```{r}
predictions <- list(prediction_1, prediction_2)
```
Now everything is ready for `enhancerDynamics()` to run.
```{r}
clusters <- crupR::getDynamics(data = predictions, C = 2)
```
The output of this step is a list containing the complete meta data file and the condition-specific clusters in a GRanges object.

##### *The parameters*

`getDynamics()` possesses a rather extensive set of paremeters that are not as self-explanatory as the ones of the previous steps. Therefore, this subsection offers more in detail explanations of their use. 

*w_0: Since comparing all the enhancer regions between two conditions to find siginificant clusters would be too expensive, a filter step was included. In this step, enhancers whose normalized prediction means did not differ strongly enough were discared before running the actual test. w_0 is the minimum difference between the means in order to be considered for the comparison. The default is 0.5.

*cutoff: The threshold value for the p-values that are computed during the clustering process. The default is 0.05.

*W: Number of bins +/- the current bin that should be included when calculating the p-values. The default is 10. However, the range for values that produce enhancer clusters of a realistic size is limited to [2, 30]. 

*C: Number of cores that should be used for parallel processing. The default is 1.

##### *Visualization of the clusters using `plotSummary()`*

As trying to decode the patterns for each cluster can be a bit complicated and also take some time, *crupR* offers the visualization function `plotSummary()` for the detected enhancer clusters. The function shows the boxplots of the median probabilities of the enhancers of each condition for each cluster. This is a simple way of seeing which clusters are active in which conditions. 

```{r}
crupR::plotSummary(clusters)
```

#### *Step 3: Find target genes of the condition-specific enhancer clusters*

The last step in the pipeline is the identification of the target genes that get regulated by the condition-specific enhancer clusters found in the previous step.


##### *Gene expression counts*
For the function `getTargets()`the output of the step before is not sufficient. *crupR* requires additionally the gene expression counts for each condition and replicate all in one GRanges object. 
*crupR* offers such an expression file for the the example sets.

```{r}
expression <- readRDS(file = system.file("extdata", "expressions.rds", package="crupR"))
```

However, there is no in-built function that calculates the gene expression counts from RNAseq experiments. That means the user has to generate this file by themselves. Please make sure that the column names for the columns that contain the gene expression counts of each sample are the same as the column names of the columns that contain the enhancer probabilities of each sample, i.e. "cond1_1", "cond1_2", "cond2_1" and so on. 

##### *Running `getTargets()`*

Looking for the enhancer targets requires - besides the gene expression counts - a .bed file containing TADs for the genome of your choice. In case mm10 is used, *crupR* provides a suitable .bed file which is also used per default for mm10 data. In this case, simply running the code below is sufficient:

```{r}
targets <- crupR::getTargets(data=clusters, expr = expression, genome = "mm10", C = 2)
```

In case another genome is used or a different BED file is preferred, providing the path to the BED file is necessary.
Per default *crupR* chooses the genes that are in the same TAD as the enhancer cluster as candidate genes. Then, the expression values of those candidate genes and enhancer probabilities are correlated to identify putative target genes. However, *crupR* also offers an alternative approach that uses the nearest gene as candidate gene. This approach does not require TAD regions and is convenient when those are not available for some reason.
The output of `getTargets()` should be a list containing two elements. The first element should be the full meta data frame that was constructed at the beginnin , because all the truncated meta data frames are merged during this step. The second element should be a GRanges file containing the computed units. There are 9 meta columns in the resulting GRanges, thus 12 columns in total:

*seqnames: chr of the dynamic enhancer region
*ranges: start and end of the dynamic enhancer region
*strand: strand of the dynamic enhancer region
*cluster: associated clusters of the dynamic enhancer region
*cond1_1 & cond1_2: the best probability values for each region per sample
*TAD_COORDINATES: the chromosome and coordinates of the TAD in which the dynamic enhancer region is located
*CORRELATED_GENE: the ID of the gene that is correlated with the dynamic enhancer region
*CORRELATED_GENE_CHR: the chromosome of the gene that is correlated with the dynamic enhancer region
*CORRELATED_GENE_PROMOTER_START: start of the promoter of the gene that is correlated with the dynamic enhancer region
*CORRELATED_GENE_PROMOTER_END: end of the promoter of the gene that is correlated with the dynamic enhancer region
*CORRELATION: correlation value


#### *Exporting the files*
After running the *crupR* pipeline, you might want to save the resulting GRanges objects. *crupR* provides the function `saveFiles()`that exports the different files in suitable formats. 

In total, there are X possible formats: "bigWig", "bedGraph", "rds", "bed", "beds" and "UCSC". However, the format to export your file to depends on the *crupR* step:

*Output of `getEnhancers()`: The GRanges file that is produced in this step can be saved as bigWig file or an .rds file. Thus, the options "bigWig" and/or "rds" can be chosen.
*Output of `getSE()`: This step produces two new GRanges files: single enhancer peak calls () and peak clusters/super enhancers (). The single enhancer peak calls can be export as a bedGraph file and the peak clusters can be exported as a BED file.
*Output of `getDynamics()`: The GRanges file produced in this step can be exported as multiple BED files ("beds"). Each bed file contains the dynamic enhancer regions of one cluster.
*Output of `getTargets()`: The GRanges file produced in this step can be exported in (UCSC) interaction format ("UCSC").

If you want to export the files, you just have to use the output list of the respetive step as input and specify your desired format(s):

```{r}
out_dir <- ""
#save the GRanges object of the getEnhancers() step
saveFiles(data = prediction_1, modes = c("bigWig", "rds"), outdir = out_dir)
#save the GRanges object of the getSE() step
saveFiles(data = se, modes = c("bedGraph", "bed"), outdir = out_dir)
#save the GRanges object of the getDynamics() step
saveFiles(data = clusters, modes = "beds", outdir = out_dir)
#save the GRanges object of the getTargets() step
saveFiles(data = targets, modes = "UCSC", outdir = out_dir)
```
