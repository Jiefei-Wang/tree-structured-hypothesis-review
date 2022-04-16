library("poolr")

BH_test <- function(pvalues){
  min(p.adjust(pvalues, "BH"))
}

stouffer_test <- function(pvalues){
  stouffer(pvalues)$p
}
tippett_test <- function(pvalues){
  tippett(pvalues)$p
}



validate_nodes <- function(taxonomy, nodes){
  root <- nodes[["root"]]
  stopifnot(length(root$leaf)==nrow(taxonomy))
  stopifnot(length(root$offspring) == length(unique(unlist(taxonomy))))
  stopifnot(length(root$children)>0)
  count <- 0
  for(i in root$children){
    count <- count + length(nodes[[i]]$leaf)
  }
  stopifnot(count==nrow(taxonomy))
}



construct_nodes <- function(taxonomy){
  if(length(unique(taxonomy[,1]))!=1)
    taxonomy <- cbind("root", taxonomy)
  nodes <- list()
  for(i in ncol(taxonomy):1){
    lvls <- unique(taxonomy[,i])
    for(lvl in lvls){
      node <- list(name = lvl)
      node$level <- i
      if(i != ncol(taxonomy)){
        node$children <- unique(taxonomy[taxonomy[,i]==lvl,i+1])
      }
      nodes <- c(nodes, list(node))
    }
  }
  nms <- sapply(nodes, function(x)x$name)
  names(nodes) <- nms
  for(i in seq_along(nodes)){
    if(!is.null(nodes[[i]]$children)){
      for(j in nodes[[i]]$children){
        nodes[[j]]$parent <- nodes[[i]]$name
      }
    }
    nodes[[i]]$leaf <- find_leaf(nodes, i)
    nodes[[i]]$offspring <- find_offspring(nodes, i)
  }
  nodes
}

find_leaf <- function(nodes, i) {
  node <- nodes[[i]]
  if(is.null(node$children)){
    return(node$name)
  }else{
    if(!is.null(node$allChildren)){
      return(node$allChildren)
    }
    children <- c()
    for(k in node$children){
      children <- c(children, find_leaf(nodes, k))
    }
    return(children)
  }
}

find_offspring <- function(nodes, i, root = TRUE) {
  node <- nodes[[i]]
  if(root)
    children <- c()
  else
    children <- node$name
  for(k in node$children){
    children <- c(children, 
                  find_offspring(nodes, k, root = FALSE)
                  )
  }
  unique(children)
}


meinshausen_FWER <- function(taxonomy, pvalues, alpha, test_fun, delta=1){
  stopifnot(length(unique(taxonomy[,1])) == 1)
  root_name <- taxonomy[1,1]
  nodes <- construct_nodes(taxonomy)
  for (i in seq_along(nodes)) {
    nodes[[i]]$pvalue <- test_fun(pvalues[nodes[[i]]$leaf])
  }
  rejections <- rep(list(NULL), ncol(taxonomy))
  names(rejections) <- colnames(taxonomy)
  if(nodes[[root_name]]$pvalue>alpha){
    return(rejections) 
  }
  rejections[[1]] <- unname(root_name)
  candidates <- nodes[[root_name]]$children
  n_leaf <- length(nodes[[root_name]]$leaf)
  while(TRUE){
    alphas <- alpha*sapply(candidates, function(x)length(nodes[[x]]$leaf))/n_leaf
    candidate_p <- sapply(candidates, function(x)nodes[[x]]$pvalue)
    idx <- which(candidate_p<alphas)
    if(length(idx)==0)
      break
    for(i in idx) {
      level <- nodes[[candidates[i]]]$level
      child_names <- nodes[[candidates[i]]]$children
      candidates <- c(candidates, child_names)
      rejections[[level]] <-c(rejections[[level]], candidates[i])
    }
    candidates <- candidates[-idx]
    if(length(candidates)==0) break
  }
  rejections
}


goeman_FWER <- function(taxonomy, pvalues, alpha, test_fun, delta=1){
  stopifnot(length(unique(taxonomy[,1])) == 1)
  root_name <- taxonomy[1,1]
  nodes <- construct_nodes(taxonomy)
  for (i in seq_along(nodes)) {
    nodes[[i]]$pvalue <- test_fun(pvalues[nodes[[i]]$leaf])
  }
  rejections <- rep(list(NULL), ncol(taxonomy))
  names(rejections) <- colnames(taxonomy)
  if(nodes[[root_name]]$pvalue>alpha){
    return(rejections) 
  }
  rejections[[1]] <- unname(root_name)
  candidates <- nodes[[root_name]]$children
  n_leaf <- length(nodes[[root_name]]$leaf)
  while(TRUE){
    alphas <- alpha*sapply(candidates, function(x)length(nodes[[x]]$leaf))/n_leaf
    candidate_p <- sapply(candidates, function(x)nodes[[x]]$pvalue)
    idx <- which(candidate_p<alphas)
    if(length(idx)==0)
      break
    for(i in idx) {
      level <- nodes[[candidates[i]]]$level
      child_names <- nodes[[candidates[i]]]$children
      candidates <- c(candidates, child_names)
      rejections[[level]] <-c(rejections[[level]], candidates[i])
      if(is.null(nodes[[candidates[i]]]$children))
        n_leaf <- n_leaf - 1
    }
    candidates <- candidates[-idx]
    if(length(candidates)==0) break
  }
  rejections
}

## use DAGmethod, but result is not consistant
goeman_FWER2 <- function(taxonomy, pvalues, alpha, test_fun){
  stopifnot(length(unique(taxonomy[,1])) == 1)
  mytest <- function(set) stouffer_test(pvalues[set])
  sets <- list()
  nms <- c()
  for(i in 1:ncol(taxonomy)){
    groups <- unique(taxonomy[,i])
    for(group in groups){
      sets <- c(sets, list(which(taxonomy[,i] == group)))
      nms <- c(nms, group)
    }
  }
  names(sets) <- nms
  struct <- suppressWarnings(construct(sets))
  
  DAG <- DAGmethod(struct, mytest, isadjusted=TRUE, alpha_max = 0.2)
  rejections <- rep(list(NULL), ncol(taxonomy))
  names(rejections) <- colnames(taxonomy)
  for(i in 1:ncol(taxonomy)){
    groups <- unique(taxonomy[,i])
    for(group in groups){
      if(!is.na(pvalue(DAG, group)) && pvalue(DAG, group)<= alpha){
        rejections[[i]] <- c(rejections[[i]], group)
      }
    }
  }
  rejections
}

Ramdas_FDR <- function(taxonomy, pvalues, alpha, test_fun) {
  nodes <- construct_nodes(taxonomy)
  root_name <- taxonomy[1,1]
  L <- nrow(taxonomy)
  p1 <- test_fun(pvalues)
  rejections <- rep(list(NULL), ncol(taxonomy))
  names(rejections) <- colnames(taxonomy)
  if(p1>alpha){
    return(rejections) 
  }
  rejections[[1]] <- unname(root_name)
  R <- c(taxonomy[1,1])
  Rc <- 1
  for(d in 2:ncol(taxonomy)){
    lvls <- unique(taxonomy[,d])
    # message("level:", d)
    
    ## find the node whose parent is significant
    included_nodes <- list()
    for(i in lvls){
      node <- nodes[[i]]
      if(node$parent%in%R){
        included_nodes <- c(included_nodes, list(node))
      }
    }
    
    included_lvls <- sapply(included_nodes, function(x)x$name)
    pvalue_level <- sapply(included_lvls, function(lvl) test_fun(pvalues[nodes[[lvl]]$leaf]))
    alpha_level <- rep(0, length(included_lvls))
    
    for(r in length(lvls):0){
      # message("progress:", 1-r/length(lvls))
      for(i in seq_along(included_nodes)) {
        node <- nodes[[i]]
        l_i <- length(node$leaf)
        m_i <- length(node$offspring) + 1
        alpha_level[i] <- alpha*l_i/L*(m_i+r+Rc-1)/m_i
      }
      if(sum(pvalue_level<=alpha_level)>=r)
        break
    }
    if(r==0)
      break
    R <- c(R, lvls[which(pvalue_level<=alpha_level)])
    Rc <- length(R)
    rejections[d] <- list(lvls[which(pvalue_level<=alpha_level)])
  }
  rejections
}

Yekutieli_FDR <- function(taxonomy, pvalues, alpha, test_fun, delta=1){
  stopifnot(length(unique(taxonomy[,1])) == 1)
  root_name <- taxonomy[1,1]
  nodes <- construct_nodes(taxonomy)
  for (i in seq_along(nodes)) {
    nodes[[i]]$pvalue <- test_fun(pvalues[nodes[[i]]$leaf])
  }
  qvalue <- alpha/2/delta
  rejections <- rep(list(NULL), ncol(taxonomy))
  names(rejections) <- colnames(taxonomy)
  if(nodes[[root_name]]$pvalue>qvalue){
    return(rejections) 
  }
  rejections[[1]] <- unname(root_name)
  R <- c(root_name)
  while(length(R)!=0){
    parent_name <- R[[1]]
    child_names <- nodes[[parent_name]]$children
    if(length(child_names)>0){
      child_pvalue <- sapply(child_names, function(x)nodes[[x]]$pvalue)
      child_pvalue_adj <- p.adjust(child_pvalue, method = "BH")
      rejected_child <- child_names[child_pvalue_adj < qvalue]
      R <- c(R[-1], rejected_child)
      lvl <- nodes[[parent_name]]$level + 1
      rejections[[lvl]] <- c(rejections[[lvl]], rejected_child)
    } else {
      R <- R[-1]
    }
  }
  rejections
}


