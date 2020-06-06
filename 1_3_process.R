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
