#本文件旨在通过前面处理好的集合，
#定义在集合的影响下的terms，proteins的语义相似性计算函数



#get_clus_ancestors将每个term的ancestors限制在每个cluster内部
#即每个cluster的所有节点与该term的ancestors做交集
#结果表示每个term在所属的cluster内部的ancestors
get_clus_ancestors <- function(clusid,
                               term,
                               m_graph = meta_graph,
                               g_graph = go_graph) {
  #cluster内部的所有节点
  clu <- m_graph[m_graph$id == clusid, ]$terms
  anc <- g_graph[g_graph$term == term, ]$ancestors

  #做交集
  intersect(unlist(anc), unlist(clu))
}




#取该term在所属的cluster(clust)下的ica值
get_com_anc_ica <- function(term, clust_id, cluster = term_cluster) {
    content <- cluster[cluster$term == term, ]
    #取对应的location
    loca <- which(unlist(content$clusid) == clust_id)
    unlist(content$ica)[loca]
}



#当term1属于clus1,term2属于clus2时，计算term1与term2的语义相似值
get_lca <- function(clus1, clus2, term1, term2) {
    if (clus1 == "meta" || clus2 == "meta") {
        return(NULL)
    }
    if (identical(clus1, clus2)) {
        #若来自同一cluster
        anc1 <- get_clus_ancestors(clus1, term1)
        anc2 <- get_clus_ancestors(clus1, term2)
        #找其在cluster内部的共同祖先
        common_anc <- intersect(anc1, anc2)
        #取对应的ica值
        value <- lapply(common_anc, get_com_anc_ica, clus1)
        value <- max(unlist(value))
  }else {
      #若来自不同cluster
      anc1 <- get_clus_ancestors("meta", clus1)
      anc2 <- get_clus_ancestors("meta", clus2)
      #找其meta分支下的共同祖先，等同原本身的共同祖先
      common_anc <- intersect(anc1, anc2)
      value <- lapply(common_anc, get_com_anc_ica, "meta")
      value <- max(unlist(value))
  }
  return(value)
}


#计算term1与term2的相似值
#因为不同ont之间的node无交集，所以不用指定ont
calculate_terms <- function(term1, term2,
                            method = "max",
                            nodes = all_node,
                            cluster = term_cluster) {
    if (!term1 %in% nodes | !term2 %in% nodes) {
        print("one of two terms has no annotations")
        return(NULL)
    } else {
        #取各自的clusters
        clus1_list <- unlist(cluster[cluster$term == term1, ]$clusid)
        clus2_list <- unlist(cluster[cluster$term == term2, ]$clusid)

        #不同的cluster下计算
        value <- sapply(clus1_list, function(e) {
            sapply(clus2_list, get_lca, e, term1, term2)
        })

      if (is.null(unlist(value))) {
          return(NULL)
      }else {
          return(operate_method(value, method))
      }
  }
}


#计算两个protein的语义相似值
protcss <- function(pro1, pro2, ont = "b",
                    method = "max",
                    pros = all_pro,
                    p_anno = pro_annotations) {
    if (!pro1 %in% pros | !pro2 %in% pros) {
        print("one of two proteins has no annotations")
        return(NULL)
    }else {
        #取各自的注释terms
        content1 <- p_anno[p_anno$ont == ont & p_anno$pro == pro1, ]
        term1_list <- unlist(content1$anno)
        content2 <- p_anno[p_anno$ont == ont & p_anno$pro == pro2, ]
        term2_list <- unlist(content2$anno)

        #两组terms遍历计算
        value <- sapply(term1_list, function(e) {
            sapply(term2_list, calculate_terms, e, method = method)
      })
    }

    if (is.null(unlist(value))) {
        NULL
    }else {
        return(operate_method(value, method))
    }
}
