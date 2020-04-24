#all the nodes will cluster into different graphs according to the cutoff

cutoff <- 3.5#ict的阈值
nodes <- names(offspring)#取所有节点，3488个

#get_log用于计算两者比值的log值
get_log <- function(item1, item2) {
  len1 <- length(item1)
  len2 <- length(item2)
  return(-log10(len1 / len2))
}

#计算ict值
ict <- lapply(offspring, get_log, nodes)
names(ict) <- nodes

#取经阈值删选过后的ict
nodes_cutoff <- names(ict[which(ict <= cutoff)])#2439个

#get_comparision用于比较两个term的ict值是否close
meta_terms <- nodes_cutoff
for (term1 in meta_terms) {
  terms <- intersect(parents[[term1]], nodes_cutoff)
  for (term2 in terms) {
    if (ict[[term2]] != 0 & ict[[term1]] / ict[[term2]] <= 1.2) {
      loca <- charmatch(term1, meta_terms)
      meta_terms <- meta_terms[-loca]
      break
    }
  }
}

#另存为data.frame格式，便于操作
meta_graph <- data.frame(id = meta_terms,
            ict = unlist(lapply(meta_terms, function(e) ict[[e]])),
            stringsAsFactors = FALSE)

#以每个term的offspring作为meta_terms的内容，组成sub-graph
#每个sub-graph去掉来自其他sub-graph的信息：后代节点
remove_dup <- function(id, meta_graph) {
  terms <- offspring[[id]]
  remove <- intersect(terms, meta_graph$id)
  remove <- setdiff(remove, id)
  remove_off <- unlist(lapply(remove, function(e) offspring[[e]]))
  terms <- setdiff(terms, remove_off)
  return(terms)
}

#meta_graph有三列，id，ict，terms
#另存为data.frame格式，便于操作
meta_graph <- data.frame(id = meta_terms,
                      ict = unlist(lapply(meta_terms, function(e) ict[[e]])),
                      terms = I(lapply(meta_graph$id, remove_dup, meta_graph)),
                      stringsAsFactors = FALSE)


#加上总的cluster meta，meta里的node即为所有的root-node
meta_graph <- rbind(meta_graph, c(id = "meta", ict = NA, terms = NA))
tmp <- list(setdiff(meta_graph$id, "meta"))
meta_graph[meta_graph$id == "meta", ]$terms <- tmp


#返回该term所属的cluster，要么是meta，要么是某个root-node
get_cluster <- function(term, meta_graph) {
    t <- lapply(meta_graph$terms, function(e) if (term %in% e) TRUE else FALSE)
    clust <- meta_graph$id[which(t == TRUE)]
    return(clust)
}

#term_cluster有两列：term，clusid（是list）
term_cluster <- data.frame(id = nodes,
                           clusid = I(lapply(nodes, get_cluster, meta_graph)),
                           stringsAsFactors = FALSE)

source("3_cluster.R")
