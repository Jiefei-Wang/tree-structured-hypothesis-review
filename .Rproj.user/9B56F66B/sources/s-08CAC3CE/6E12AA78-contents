library(cherry)
library("BOUTH")
source("read_data.R")
source("functions.R")

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


a=test.1$results.by.node[test.1$results.by.node$is.detected, ]
a
a[which(a$otu==" "),]$family
result_Yekutieli$family[!result_Yekutieli$family%in%a[which(a$otu==" "),]$family]

select <- "Firmicutes_Clostridia_[Mogibacteriaceae]"
result_Yekutieli$otu[startsWith(result_Yekutieli$otu, select)]
a$otu[startsWith(a$otu, select)]
dim(taxonomy[taxonomy$family==select,])





library("igraph")
library("ape")
library(structSSI)
nodes <- construct_nodes(taxonomy)
edges <- c()
tree.p <- c()
for(i in rev(seq_along(nodes))){
  node <- nodes[[i]]
  tree.p <- c(tree.p, mytest(node$leaf))
  for(j in node$children){
    edges <- c(edges,node$name, j)
  }
}
names(tree.p)<-rev(names(nodes))


mygraph <- make_graph(edges)
tree.el <- get.edgelist(mygraph)
tree.p <- tree.p[V(mygraph)$name]
hyp.tree <- hFDR.adjust(tree.p, tree.el, 0.2)
plot(hyp.tree, adjust = FALSE)

summary(hyp.tree)


# We can visualize the difference between the unadjusted and the adjusted
# trees.
plot(hyp.tree, adjust = FALSE)
plot(hyp.tree, adjust = TRUE)



# 
# data(IBD)
# dim(IBD$tax.table)
# colnames(IBD$tax.table)
# filter <- which(IBD$tax.table$phylum=="Firmicutes")
# taxonomy <- IBD$tax.table[filter,c("phylum", "species", "otu")]
# pvalues <- IBD$pvalue.otus[filter]
# taxonomy1 <- cbind("root", taxonomy)
# IBD$tax.table[IBD$tax.table$species=="durum",]
# taxonomy <- IBD$tax.table
# pvalues <- IBD$pvalue.otus
