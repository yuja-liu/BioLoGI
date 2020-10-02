setwd("~/YABMSS")
library(DrImpute)
source("scripts/preprocess.R")    # Curiously, func preprocess cannot be imported

# Load read count of GSE75748
tpm <- read.csv("data/GSE75748/hg38_GSE75748.tpm", sep = '\t', header = TRUE)
head(tpm[,1:10])

# Imputation
tpm_h9.mat = data.matrix(tpm[, 3:164])    # select only h9 section
rownames(tpm_h9.mat) <- tpm$id
tpm_h9.mat <- preprocess(tpm_h9.mat, min.expressed.gene = 0)    # expressed in only 0 or 1 cell is removed
print(paste(nrow(tpm_h9.mat), "genes remains"))
tpm_h9.log <- log(tpm_h9.mat + 1)
tpm_h9.imp <- DrImpute(tpm_h9.log, mc.cores = 4)
head(tpm_h9.imp[, 1:10])

# save to textfile
write.table(tpm_h9.imp, file = "data/GSE75748/h9_imputed.tsv", quote = FALSE, sep = '\t')
