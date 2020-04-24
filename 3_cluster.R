#term_annotations是每个term对应的genes信息
term_annotations <- final_annotations

#meta_graph以每个root-node为key
#get_max_ica计算每个root里面最大的ica值
#ica = 该term的genes个数 / 总genes数目
#总genes数目等于BP("GO:0008150")本身对应的genes数目 ：6357
get_max_ica <- function(clusid, meta_graph, term_annotations) {
  all_node <- unlist(meta_graph[meta_graph$id == clusid, ]$terms)
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
clus_id <- meta_graph$id
meta_graph$max_ica <- lapply(clus_id, get_max_ica, meta_graph, term_annotations)

#get_clus_ica计算每个term在不同的cluster里不同的ica值
#等于 该term本身的ica / 所属cluster里的max_ica
get_clus_ica <- function(term, term_annotations, meta_graph) {
  clusid_list <- unlist(term_cluster[term_cluster$id == term, ]$clusid)
  #该term所属的所有clusters
  own_ica <- get_log(term_annotations[[term]], term_annotations[["GO:0008150"]])
  #本身原本的ica值
  clus_max_ica <- lapply(clusid_list,
                         function(e) meta_graph[meta_graph$id == e, ]$max_ica)
  #所有cluster的max_ica
  clus_ica <- lapply(unlist(clus_max_ica), function(e) own_ica / e)
  #本身的ica除以每个cluster的max_ica
  return(unlist(clus_ica))
}

#term_cluster共三列：term, clusid, ica
term_cluster$ica <- lapply(term_cluster$id, get_clus_ica,
                           term_annotations,
                           meta_graph)

#对每一个gene有其对应的term列表，其中应去掉内部的父子关系(ancestor)
gene_remove_redun <- function(term, term_list, offs = offspring) {
  term_list <- setdiff(term_list, term)
  result <- intersect(offs[[term]], term_list)
  if (length(result) == 0) {
    return(term)
  }else {
    return(NULL)
  }
}

#gene_annnotations
genes <- unique(gene_term$gene)

#从变量gene_term中提取出每一个gene的注释信息（term列表）
get_gene_anno <- function(gene, gene_terms = gene_term) {
  term_list <- gene_terms[gene_terms$gene == gene, ]$term
  term_list <- intersect(term_list, nodes)
  term_list <- unlist(lapply(term_list, gene_remove_redun, term_list))
  return(term_list)
}

gene_annotations <- lapply(genes, get_gene_anno)
names(gene_annotations) <- genes

source("4_calculate.R")
