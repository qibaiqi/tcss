#用以根据以上信息处理每个蛋白的注释信息
#即每个蛋白对应的term_list，需要去冗余




#pro_annnotations
#所有的proteins
all_pro <- unique(pro_term$pro)

#对term列表去冗余：去掉内部的父子关系(offspring)
#判断该term在term_list中是否有父子关系
pro_remove_redun <- function(term, ont, term_list, graph = go_graph) {
  #先去掉其本身
  term_list <- setdiff(term_list, term)
  #该term的后代节点
  offs <- graph[graph$ont == ont & graph$term == term, ]$offspring
  #做交集
  result <- intersect(unlist(offs), term_list)
  #为空则表明没有，可留下该term
  if (length(result) == 0) {
    return(term)
  }
}



#从变量pro_term中提取出每一个pro的注释信息:term列表
#从不同的ontology提取到的结果不同，所以表明ont
get_pro_anno <- function(pro, ont, g_graph = go_graph, pro_terms = pro_term) {
  #从pro_term中提取
  terms <- pro_term[pro_term$pro == pro1, ][, "term"]
  #与该ont的节点做交集
  term_list <- intersect(terms, g_graph[g_graph$ont == ont, ][, "term"])
  #去冗余
  unlist(lapply(term_list, pro_remove_redun, ont, term_list))
}



#共三列:protein, ontology, annotations
pro_annotations <- data.frame(pro = rep(all_pro, each = 3),
                              ont = rep(c("b", "m", "c"), length(all_pro)),
                              stringsAsFactors = F)

pro_annotations$anno <- mapply(get_pro_anno,
                               pro_annotations$pro,
                               pro_annotations$ont)
