#综合两份源文件，只留下注释为空的节点及其所有信息
#得到 ancestors, offspring, 蛋白注释信息: pro_annotations



#depth_searth是迭代函数，找到parent的parent, children的children
depth_search <- function(term, pool) {
    items <- pool[[term]]#该term的parents或children
    new_items <- unlist(lapply(items, function(e) pool[[e]]))
    items <- unique(c(new_items, items))#合并并去重
    if (!is.null(new_items)) {
        new_items <- unlist(lapply(new_items, depth_search, pool))
        items <- unique(c(new_items, items))
    }
    items <- c(items, term)#包括其本身
}



#children的所有children...迭代得到 offspring
offspring <- lapply(total_node$id, depth_search, children)
names(offspring) <- total_node$id


#annotations信息的合并

#父节点的注释信息等于
#其本身的注释信息及其所有后代节点的注释信息的总和
anno_add <- function(term, offs = offspring, annotations = init_annotations) {
    off_list <- offs[[term]]
    anno_list <- lapply(off_list, function(e) annotations[[e]])
    anno_list <- c(unlist(anno_list), annotations[[term]])
    unique(anno_list)
}

#最终的annotations
final_annotations <- lapply(total_node$id, anno_add)
names(final_annotations) <- total_node$id


#各ontology的根节点, BP, MF, CC
root_node <- c("GO:0008150", "GO:0003674", "GO:0005575")

#作者认为根节点的所有后代节点为ont的所有节点
#三个ont的放在一起
ont_node <- lapply(root_node, depth_search, children)


#从ont_node中留下final_annotations不为空的节点
node_with_anno <- which(lapply(final_annotations, length) != 0)
ont_node <- lapply(ont_node, intersect, total_node$id[node_with_anno])
names(ont_node) <- c("b", "m", "c")



all_node <- unlist(ont_node)

#迭代得到ancestors, offspring后，在其内部信息也要过滤一遍
get_anc_off <- function(term, pool, nodes = all_node) {
    content <- depth_search(term, pool)
    #与有注释的节点做交集
    intersect(content, nodes)
}


#go_graph包含所有term的ont, ancestors, offspring, annotations信息
go_graph <- data.frame(term = all_node,
                       ont = unlist(lapply(all_node, function(e) {
                           total_node[total_node$id == e, ]$ont
                       })),
                       ancestors = I(lapply(all_node, get_anc_off, parents)),
                       offspring = I(lapply(all_node, get_anc_off, children)),
                       annotations = I(lapply(all_node, function(e) {
                           final_annotations[[e]]
                       })),
                       stringsAsFactors = F)








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
get_pro_anno <- function(pro, ont, graph = go_graph, pro_terms = pro_term) {
    #从pro_term中提取
    term_list <- pro_terms[pro_terms$pro == pro, ]$term
    #与该ont的节点做交集
    term_list <- intersect(term_list, graph[graph$ont == ont, ]$term)
    #去冗余
    term_list <- unlist(lapply(term_list, pro_remove_redun, ont, term_list))
}


#共三列:protein, ontology, annotations
pro_annotations <- data.frame(pro = rep(all_pro, each = 3),
                              ont = rep(c("b", "m", "c"), length(all_pro)),
                              stringsAsFactors = F)

pro_annotations$anno <- mapply(get_pro_anno,
                              pro_annotations$pro,
                              pro_annotations$ont)
source("2_1_clustering.R")