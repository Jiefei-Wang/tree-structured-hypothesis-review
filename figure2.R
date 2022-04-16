library(BOUTH)
source("functions.R")

## Create the taxonomy tree
taxonomy <- 
  cbind(
    rep("H11", 5),
    c("H21","H21","H22","H22","H22"),
    paste0("H3",1:5)
  )
colnames(taxonomy) <- paste0("level", 1:3)

## pvalues for the leaf hypotheses
pvalues <- c(0.03,0.03,0.01,0.03,0.03)
names(pvalues) <- paste0("H3",1:5)

## check the p-value of the root hypothesis
## Not match!!!
tippett_test(pvalues)

## A dummy global test for reproducing the result in Figure 1
dummy_test <- function(pvalues){
  if(length(pvalues) == 1)
    return(pvalues)
  if(all(names(pvalues) %in% c("H31","H32")))
    return(0.045)
  if(all(names(pvalues) %in% c("H33","H34","H35")))
    return(0.055)
  return(0.02)
}


## Choose a global test
test_func <- dummy_test
alpha <- 0.1

##FWER
result_meinshausen <- meinshausen_FWER(taxonomy, pvalues, alpha=alpha, test_fun = test_func)
result_meinshausen

result_goeman <- goeman_FWER(taxonomy, pvalues, alpha=alpha, test_fun = test_func)
result_goeman

# FDR
result_Ramdas <- Ramdas_FDR(taxonomy, pvalues, alpha=alpha, test_fun = test_func)
result_Ramdas

result_Yekutieli <- Yekutieli_FDR(taxonomy, pvalues, alpha=alpha, test_fun = test_func)
result_Yekutieli

# FAR
test.1 <- bouth(anno.table = taxonomy, pvalue.leaves = pvalues,
                na.symbol = "unknown", far = 0.1, is.weighted = TRUE)

test.1$results.by.level

