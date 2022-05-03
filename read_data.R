library("BOUTH")
data(IBD)

# filter <- which(IBD$tax.table$phylum=="Firmicutes")
# filter <- 1:nrow(IBD$tax.table)
# taxonomy <- IBD$tax.table[filter,c("phylum", "class","family", "otu")]
# taxonomy <- IBD$tax.table[filter,c("phylum", "class","family", "otu")]
# pvalues <- IBD$pvalue.otus[filter]


taxonomy <- IBD$tax.table
pvalues <- IBD$pvalue.otus

taxonomy <- taxonomy[, colnames(taxonomy)!="species"]
colnames(taxonomy)[colnames(taxonomy) == "otu"] <- "species"

for(i in 2:(ncol(taxonomy))){
  taxonomy[,i] <- paste0(taxonomy[,i-1],"_",taxonomy[,i])
}
names(pvalues) <- taxonomy$species




# for(i in 1:(length(sets)-1)){
#   current <- sets[[i]]
#   for(j in (i+1):length(sets)){
#     current2<-sets[[j]]
#     stopifnot(!setequal(current,current2))
#   }
# }



