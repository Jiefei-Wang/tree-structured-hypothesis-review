library("BOUTH")
data(IBD)

filter <- which(IBD$tax.table$phylum=="Firmicutes")
taxonomy <- IBD$tax.table[filter,c("phylum", "class","family", "otu")]
pvalues <- IBD$pvalue.otus[filter]

for(i in 2:(ncol(taxonomy))){
  taxonomy[,i] <- paste0(taxonomy[,i-1],"_",taxonomy[,i])
}
names(pvalues) <- taxonomy$otu




# for(i in 1:(length(sets)-1)){
#   current <- sets[[i]]
#   for(j in (i+1):length(sets)){
#     current2<-sets[[j]]
#     stopifnot(!setequal(current,current2))
#   }
# }



