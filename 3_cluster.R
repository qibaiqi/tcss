#term_annotations是每个term对应的genes信息
term_annotations <- final_annotations

#select函数用以从data.frame中取数据
select <- function(target, from, where, wanted) {
  result <- from[from[where] == target, ][wanted]
  result <- unique(unlist(result))
  return(result)
}

#meta_terms以每个root-node为key
#get_max_ica计算每个root里面最大的ica值
#ica = 该term的genes个数 / 总genes数目
#总genes数目等于BP("GO:0008150")本身对应的genes数目 ：6357
get_max_ica <- function(clusid) {
  all_node <- select(clusid, meta_terms, "id", "terms")
  #该cluster里的所有节点
  all_node_anno <- lapply(all_node, function(e) term_annotations[[e]])
  #取所有节点的genes信息
  ica_value <- lapply(all_node_anno, get_log, term_annotations[["GO:0008150"]])
  #分别计算所有节点的ica
  max_ica <- max(unlist(ica_value))
  #选取最大值
  return(max_ica)
}

#meta_terms有四列：id, ict, terms, max_ica
meta_terms$max_ica <- lapply(meta_terms$id, get_max_ica)

#get_clus_ica计算每个term在不同的cluster里不同的ica值
#等于 该term本身的ica / 所属cluster里的max_ica
get_clus_ica <- function(term, clusid_list) {
  clusid_list <- unlist(clusid_list)
  #该term所属的所有clusters
  own_ica <- get_log(term_annotations[[term]], term_annotations[["GO:0008150"]])
  #本身原本的ica值
  clus_max_ica <- lapply(clusid_list, select, meta_terms, "id", "max_ica")
  #所有cluster的max_ica
  clus_ica <- lapply(unlist(clus_max_ica), function(e) own_ica / e)
  #本身的ica除以每个cluster的max_ica
  return(unlist(clus_ica))
}

#term_cluster共三列：term, clusid, ica
term_cluster$ica <- mapply(get_clus_ica, term_cluster$term, term_cluster$clusid)

#get_clus_ancestors将每个term的ancestors限制在每个cluster内部
#即每个cluster的所有节点与该term的ancestors做交集
#结果表示每个term在所属的cluster内部的ancestors
get_clus_ancestors <- function(clusid, term) {
  clu <- meta_terms[meta_terms$id == clusid, ]$terms
  #该cluster内的所有节点
  clu_anc <- intersect(ancestors[[term]], unlist(clu))
  #做交集
  return(unlist(clu_anc))
}

#最终形成graph变量,是一个列表
#graph:
#[[term]]...
#       [[clusid]]
#                [[ancestors]]...
#       [[ica]]...
#       [[ancestors]]...
#       [[genes]]...

graph <- list()

for (i in seq_len(nrow(term_cluster))) {
  term <- term_cluster$term[i]
  #第一层
  clusid_list <- term_cluster$clusid[i]
  #该term对应的clusters
  graph[[term]][["clusid"]] <- clusid_list
  #第二层
  tmp <- lapply(unlist(clusid_list), get_clus_ancestors, term)
  #cluster内部的ancestors
  graph[[term]][["clusid"]][["ancestor"]] <- tmp
  #第三层
  graph[[term]][["ica"]] <- term_cluster$ica[i]
  #第二层
  graph[[term]][["ancestors"]] <- ancestors[[term]]
  #第二层
  graph[[term]][["genes"]] <- term_annotations[[term]]
  #第二层
}
