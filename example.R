library(cherry)
library("BOUTH")
source("read_data.R")
source("functions.R")

lengths(apply(taxonomy, 2, function(x) unique(x)))

#FWER
result_meinshausen <- meinshausen_FWER(taxonomy, pvalues, alpha=0.2, test_fun = stouffer_test)
lengths(result_meinshausen)

result_goeman <- goeman_FWER(taxonomy, pvalues, alpha=0.2, test_fun = stouffer_test)
lengths(result_goeman)

# FDR
result_Ramdas <- Ramdas_FDR(taxonomy, pvalues, alpha=0.2, test_fun = stouffer_test)
lengths(result_Ramdas)

result_Yekutieli <- Yekutieli_FDR(taxonomy, pvalues, alpha=0.2, test_fun = stouffer_test)
lengths(result_Yekutieli)

# FAR
test.1 <- bouth(anno.table = taxonomy, pvalue.leaves = pvalues,
               na.symbol = "unknown", far = 0.1, is.weighted = TRUE)

test.1$results.by.level

