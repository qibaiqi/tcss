#用以求出各个sub-root-node所对应的cluster内部的最大的ica值
#再据此求出各个term对应的修正后的ica值


#计算该term的ica值, information content value, 作为分子
get_ica <- function(term, g_graph = go_graph, root_nodes = root_node) {
    content1 <- g_graph[g_graph$term == term, ]
    #该term的ont
    ont <- content1$ont
    #该term的对应蛋白数目
    len1 <- length(unlist(content1[["annotations"]]))
    #取各个ont的根节点的蛋白信息
    content2 <- g_graph[charmatch(root_nodes, g_graph$term), ][["annotations"]]
    #取各个ont的的的蛋白数目
    len2 <- switch(ont,
                   "b" = length(unlist(content2[1])),
                   "m" = length(unlist(content2[2])),
                   "c" = length(unlist(content2[3])))

  - log10(len1 / len2)
}



#取该cluster内部最大的ica值，作为分母
#meta_graph以每个root-node为key
#ica = 该term的pros个数 / 总pros数目
get_max_ica <- function(clusid, ont,
                        m_graph = meta_graph,
                        annotations = term_annotations) {
    content <- m_graph[m_graph$id == clusid & m_graph$ont == ont, ]
    #该cluster内部的所有节点
    all_node <- unlist(content$terms)
    #取其ica值
    ica_value <- sapply(all_node, get_ica)
    #取最大
    max(ica_value)
}



#第4列
meta_graph$max_ica <- mapply(get_max_ica, meta_graph$id, meta_graph$ont)





#get_clus_ica计算每个term在不同的cluster里不同的ica值
#修正ica = 该term本身的ica / 所属cluster里的max_ica
get_clus_ica <- function(term, m_graph = meta_graph, cluster = term_cluster) {
    content1 <- cluster[cluster$term == term, ]
    #该term的ont, 所属clusters, 原本ica值
    ont <- content1$ont
    clus_list <- unlist(content1$clusid)
    own_ica <- get_ica(term)

    #所属clusters的各自max_ica
    content2 <- m_graph[m_graph$ont == ont, ]
    clus_max_icas <- content2[charmatch(clus_list, content2$id), ]$max_ica
    #修正ica
    clus_ica <- unlist(lapply(clus_max_icas, function(e)
        if (e != 0) {
            own_ica / e
        } else 0))
    return(clus_ica)
}





#返回该term所属的cluster，要么是meta，要么是某个sub-root-node
get_cluster <- function(term, graph = meta_graph) {
    result <- lapply(graph$terms, function(e) term %in% e)
    graph$id[which(result == TRUE)]
}


#term_cluster有4列：term, ont, clusid(list), ica(list)
term_cluster <- go_graph[c("term", "ont")]

term_cluster$clusid <- lapply(term_cluster$term, get_cluster)

term_cluster$ica <- lapply(term_cluster$term, get_clus_ica)
