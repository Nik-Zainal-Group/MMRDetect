# MMRDetect

## Overview

MMRDetect is a mutational signature-based classifier for identifying tumors with mismatch repair deficiency using whole-genome sequencing data.


## Installation

```{r, eval = FALSE}
# Install the released version from Github
git clone https://github.com/xqzou/MMRDetect.git
cd MMRDetect
R CMD INSTALL .
```


## Example
Apply MMRDetect to 26 breast cancers

```{r, eval = FALSE}
# import the package 
library(MMRDetect)

# 1) Generate substitution catalogue
all_subs <- read.table("./breast26_subs.txt",sep = "\t", header = T, as.is = T)
sub_catalogues <- GenCatalogue(all_subs,"Sample")

# 2) Generate indel catalogue
all_indels <- read.table("./breast26_indels.txt",sep = "\t", header = T, as.is = T)
indel_classied <- indel_classifier(all_indels)
indel_catalogues <- gen_indelmuttype_MMRD(indel_classied, "Sample","indeltype_short")

# 3) Calculate variables for MMRDetect
muts_variables <- MMRDetect.compute.variables(sub_catalogues, indel_catalogues, "Breast")
```
You can find output file of calculated variables `sample_feature_summary_Breast.txt`.

```{r, eval = FALSE}

# 4) Apply MMRDetect
MMRDetect_classified <- MMRDetect.classify(muts_variables)
write.table(MMRDetect_classified,"MMRDetect_classified.txt",sep="\t", col.names = T, row.names = F, quote = F)
```


