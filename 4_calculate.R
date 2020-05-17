#本文件旨在通过前面处理好的集合，
#定义在集合的影响下的terms，proteins的语义相似性计算函数

#get_clus_ancestors将每个term的ancestors限制在每个cluster内部
#即每个cluster的所有节点与该term的ancestors做交集
#结果表示每个term在所属的cluster内部的ancestors
get_clus_ancestors <- function(clusid,
                               term,
                               graph = meta_graph,
                               anc = ancestors) {
    clu <- graph[graph$id == clusid, ]$terms
    clu_anc <- intersect(anc[[term]], unlist(clu))
    return(unlist(clu_anc))
}

#取该term在所属的cluster(clust)下的ica值
get_com_anc_ica <- function(term, clust, cluster = term_cluster) {
    content <- subset(cluster, id == term)
    loca <- which(unlist(content$clusid) == clust)
    value <- unlist(content$ica)[loca]
}


#当term1属于clus1,term2属于clus2时，计算term1与term2的语义相似值
get_lca <- function(clus1, clus2, term1, term2) {
    if (clus1 == "meta" || clus2 == "meta") {
      return(NULL)
    }
    if (identical(clus1, clus2)) {
        anc1 <- get_clus_ancestors(clus1, term1)
        anc2 <- get_clus_ancestors(clus1, term2)
        common_anc <- intersect(anc1, anc2)
        value <- unlist(lapply(common_anc, get_com_anc_ica, clus1))
    }else {
        anc1 <- get_clus_ancestors("meta", clus1)
        anc2 <- get_clus_ancestors("meta", clus2)
        common_anc <- intersect(anc1, anc2)
        value <- unlist(lapply(common_anc, get_com_anc_ica, "meta"))
    }
    return(value)
}


#计算term1与term2的相似值，取最大值
calculate_terms <- function(term1, term2,
                            all_nodes = nodes,
                            cluster = term_cluster) {
    if (!term1 %in% all_nodes | !term2 %in% all_nodes) {
        print("one of two terms has no annotations")
        return(NULL)
    }else {
        clus1_list <- unlist(subset(cluster, id == term1)$clusid)
        clus2_list <- unlist(subset(cluster, id == term2)$clusid)
        value <- unlist(lapply(clus1_list, function(e){
            lapply(clus2_list, get_lca, e, term1, term2)
        }))
        if (is.null(value)) {
            return(NULL)
        }else {
            return(max(value))
            }
        }
}


#计算两个protein的语义相似值
protcss <- function(pro1, pro2,
                    all_pros = pros,
                    pro_anno = pro_annotations) {
    if (!pro1 %in% all_pros | !pro2 %in% all_pros) {
        print("one of two proteins has no annotations")
        return(NULL)
    }else {
        term1_list <- pro_anno[[pro1]]
        term2_list <- pro_anno[[pro2]]
        value <- unlist(lapply(term1_list, function(e){
            lapply(term2_list, calculate_terms, e)
        }))
        if (is.null(value)) {
            return(NULL)
        }else {
            return(max(value))
            }
        }
}

source("5_roc.R")
