#all the nodes will cluster into different graphs according to the cutoff
#本文件旨在根据每个term的ict值将所有的term分类，进而根据父子关系形成不同的集合

cutoff <- 4.0#ict的阈值
nodes <- names(offspring)#取所有节点，3488个

#get_log用于计算两者比值的log值
get_log <- function(item1, item2) {
  -log10(length(item1) / length(item2))
}

#计算ict值
ict <- lapply(offspring, get_log, nodes)
names(ict) <- nodes

#取经阈值删选过后的ict
nodes_cutoff <- names(ict[which(ict <= cutoff)])#2439个

#get_comparision用于比较两个term的ict值是否close,是则删除
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
remove_dup <- function(id, offs = offspring, graph = meta_graph) {
    terms <- offs[[id]]
    remove <- intersect(terms, graph$id)
    remove <- setdiff(remove, id)
    remove_off <- unlist(lapply(remove, function(e) offs[[e]]))
    terms <- setdiff(terms, remove_off)
    return(terms)
}

#meta_graph有三列，id，ict，terms
#另存为data.frame格式，便于操作
meta_graph <- data.frame(id = meta_terms,
                      ict = unlist(lapply(meta_terms, function(e) ict[[e]])),
                      terms = I(lapply(meta_graph$id, remove_dup)),
                      stringsAsFactors = FALSE)


#加上总的cluster:meta，meta里的node即为所有的root-node
meta_graph <- rbind(meta_graph, c(id = "meta", ict = NA, terms = NA))
tmp <- list(setdiff(meta_graph$id, "meta"))
meta_graph[meta_graph$id == "meta", ]$terms <- tmp


#返回该term所属的cluster，要么是meta，要么是某个root-node
get_cluster <- function(term, graph = meta_graph) {
    result <- lapply(graph$terms, function(e) if (term %in% e) TRUE else FALSE)
    graph$id[which(result == TRUE)]
}

#term_cluster有两列：term，clusid（是list）
term_cluster <- data.frame(id = nodes,
                           clusid = I(lapply(nodes, get_cluster)),
                           stringsAsFactors = FALSE)

source("3_cluster.R")
