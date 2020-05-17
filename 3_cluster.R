#本文件旨在计算统计并计算每个集合的root-node，其ict值，其所含的所有terms，
#及集合内部最大的ica值。同时计算每个term所属集合，及不同的集合下对应的ica值
#并且根据前面源文件处理形成每个protein的注释列表信息

#term_annotations是每个term对应的proteins信息
term_annotations <- final_annotations

#meta_graph以每个root-node为key
#get_max_ica计算每个root里面最大的ica值
#ica = 该term的pros个数 / 总pros数目
#总pros数目等于BP("GO:0008150")本身对应的pros数目 ：6357
get_max_ica <- function(clusid,
                        graph = meta_graph,
                        annotations = term_annotations) {
    all_node <- unlist(graph[graph$id == clusid, ]$terms)
    #该cluster里的所有节点
    all_node_anno <- lapply(all_node, function(e) annotations[[e]])
    #取所有节点的pros信息
    ica_value <- lapply(all_node_anno, get_log, annotations[["GO:0008150"]])
    #分别计算所有节点的ica,GO:0008150是BP节点本身，他的注释是所有的proteins
    max_ica <- max(unlist(ica_value))
    #选取最大值
}

#meta_terms有四列：id, ict, terms, max_ica
clus_id <- meta_graph$id
meta_graph$max_ica <- lapply(clus_id, get_max_ica)

#get_clus_ica计算每个term在不同的cluster里不同的ica值
#等于 该term本身的ica / 所属cluster里的max_ica
get_clus_ica <- function(term,
                         annotations = term_annotations,
                         graph = meta_graph,
                         cluster = term_cluster) {
    clusid_list <- unlist(cluster[cluster$id == term, ]$clusid)
    #该term所属的所有clusters
    own_ica <- get_log(annotations[[term]], annotations[["GO:0008150"]])
    #本身原本的ica值
    clus_max_ica <- lapply(clusid_list,function(e){
        graph[graph$id == e, ]$max_ica
    })
    #所有cluster的max_ica
    clus_ica <- lapply(unlist(clus_max_ica), function(e) own_ica / e)
    #本身的ica除以每个cluster的max_ica
    return(unlist(clus_ica))
}

#term_cluster共三列：term, clusid, ica
term_cluster$ica <- lapply(term_cluster$id, get_clus_ica)

#对每一个protein有其对应的term列表，其中应去掉内部的父子关系(ancestor)
pro_remove_redun <- function(term, term_list, offs = offspring) {
    term_list <- setdiff(term_list, term)
    result <- intersect(offs[[term]], term_list)
    if (length(result) == 0) {
        return(term)
    }else {
        return(NULL)
        }
}

#pro_annnotations
#所有的proteins
pros <- unique(pro_term$pro)

#从变量pro_term中提取出每一个pro的注释信息（term列表）
get_pro_anno <- function(pro, all_nodes = nodes, pro_terms = pro_term) {
    term_list <- pro_terms[pro_terms$pro == pro, ]$term
    term_list <- intersect(term_list, all_nodes)
    term_list <- unlist(lapply(term_list, pro_remove_redun, term_list))
}

pro_annotations <- lapply(pros, get_pro_anno)
names(pro_annotations) <- pros

source("4_calculate.R")
