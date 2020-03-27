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
ict <- lapply(children, get_log, children)
names(ict) <- nodes

#取经阈值删选过后的ict
nodes_cutoff <- ict[which(ict <= cutoff)]

#另存为data.frame格式，便于操作
meta_terms <- data.frame(id = names(nodes_cutoff),
                         ict = unlist(unname(nodes_cutoff)),
                         stringsAsFactors = FALSE)

#ciduanwuyong？？？
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
#但是就这样直接删除，不需要合并吗？不然一些node没有了对应的root节点
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

#以每个term的children作为meta_terms的内容
meta_terms$terms <- lapply(meta_terms$id, function(e) children[[e]])


#if edges a->b->c and a->c exist then a->c is removed.
#此段没有运行!!! 还不是很确定
#get_relation中 term1是 a, all_node 是 b,c,d,,,判断是否存在上述情况
#get_relation <- function(term1, all_node) {
#  unlist(lapply(all_node, function(term2) {
#  if (term1 != term2 & term2 %in% offspring[[term1]]) {
#    #确定是offspring吗？
#    return(term2)
#  }else{
#    return(NULL)
#  }
#  }))
#}



#remove_unwanted 对所有的meta_terms一个个检查，看“他们之间”是否存在上述情况
#若有，需从parents,children(??)中删去那些meta_terms
#remove_unwanted <- function(item, all_node) {
#  child <- children[[item]]
#  child_meta <- intersect(child, all_node)
#  remove_term <- unlist(lapply(child_meta, get_relation, child_meta))
#  if (!is.null(remove_term)) {
#    return(remove_term)
#    children[[item]] <- setdiff(children[[item]], remove_term)
#    for (rt in remove_term) {
#      parents[[rt]] <- setdiff(parents[[rt]], item)
#    }
#  } 
#}

#lapply(meta_terms$id, remove_unwanted, meta_terms$id)

#clustering and calculate ICA
#go_annotations
#[cluster]:[clusid][ancestors][ICA]

get_cluster <- function(term, meta_terms) {
  if (term %in% meta_terms$id) {
    return(c("meta"))
  }else {
    t <- lapply(meta_terms$terms, function(e) if (term %in% e) TRUE else FALSE)
    clust <- meta_terms$id[which(t == TRUE)]
    return(clust)
  }
}

term_cluster <- data.frame(term = nodes)
term_cluster$clusid <- lapply(nodes, get_cluster, meta_terms)
