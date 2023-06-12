# **crupR - Short Description**   

crupR is improved version of the enhancer detection pipeline CRUP ([C]ondition-specific [R]egulatory [U]nits [P]rediction). Its workflow consists of the same three main steps as the original pipeline (enhancer prediction, condition-specific enhancer detection, target gene detection) and an additional pre-preapring step (normalization). It uses different layers of epigenetic information in the form of ChIP-seq data from histone modification to provide a comprehensive list of regulatory units
consisting of dynamically changing enhancers and target genes.

### *Installation*
Install crupR using following commad:

```
devtools::install_git("https://github.com/akbariomgba/crupR")
```


#### *Contact*

omgba@molgen.mpg.de

#### *Citation*

Since the crupR application note has not been submitted yet, you can cite this repository and the original paper if you use crupR:

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1860-7
