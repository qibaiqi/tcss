#all the nodes will cluster into different graphs according to the cutoff

cutoff <- 3.5#ict的阈值
nodes <- names(offspring)#取所有节点，3488个

#get_log用于计算两者比值的log值
get_log <- function(item1, item2) {
  len1 <- length(item1)
  len2 <- length(item2)
  return(-log10(len1/len2))
}

#计算ict值
ict <- lapply(offspring, get_log, offspring)
names(ict) <- nodes

#取经阈值删选过后的ict
nodes_cutoff <- ict[which(ict <= cutoff)]

#另存为data.frame格式，便于操作
meta_terms <- data.frame(id = names(nodes_cutoff),
                         ict = unlist(unname(nodes_cutoff)),
                         stringsAsFactors = FALSE)

#get_comparision用于比较两个term的ict值是否close
get_comparision <- function(term2, term1) {
  ict1 <- ict[[term1]]
  ict2 <- ict[[term2]]
  if (ict2 != 0 & ict1/ict2 < 1.2) {
    return(TRUE)
  }else {
    return(FALSE)
  }
}

#If the ict values of parent - child terms is
#in close proximity of each other then the child term is removed.
remove_close_terms <- function(term, all_nodes) {
  terms2 <- intersect(parents[[term]], all_nodes)
  if (!is.null(terms2)) {
    result <- lapply(terms2, get_comparision, term)
    if ("TRUE" %in% result) {
      return(NULL)
    }else{
      return(term)
    }
  }
}

#从meta_terms中删除那些判定为close的term
remove_terms <- lapply(meta_terms$id, remove_close_terms, meta_terms$id)
meta_terms <- meta_terms[-which(unlist(remove_terms) %in% meta_terms$id), ]

#以每个term的offspring作为meta_terms的内容
meta_terms$terms <- lapply(meta_terms$id, function(e) offspring[[e]])
