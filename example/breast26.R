# Apply MMRDetect to 26 breast cancers

library(MMRDetect)
# 1) Generate substitution catalogue
all_subs <- read.table("./breast26_subs.txt",sep = "\t", header = T, as.is = T)
sub_catalogues <- GenCatalogue(all_subs,"Sample")

# 2) Generate indel catalogue
all_indels <- read.table("./breast26_indels.txt",sep = "\t", header = T, as.is = T)
indel_classied <- indel_classifier(all_indels)
indel_catalogues <- gen_indelmuttype_MMRD(indel_classied, "Sample","indeltype_short")

# 3) Generate variables
muts_variables <- MMRDetect.compute.variables(sub_catalogues, indel_catalogues, "Breast")

# 4) Apply MMRDetect
MMRDetect_classified <- MMRDetect.classify(muts_variables)
