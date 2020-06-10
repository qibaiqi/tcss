#all the nodes will cluster into different graphs according to the cutoff
#本文件旨在根据每个term的ict值将所有的term分类
#进而根据父子关系形成不同的集合



#计算ICT的函数，ICT只用于分割为clusters
get_ict <- function(term, g_graph = go_graph, nodes = ont_node) {
    content <- g_graph[g_graph$term == term, ]
    #该节点后代数目
    len1 <- length(unlist(content$offspring))
    #该ont里所有节点数目
    len2 <- length(nodes[[content$ont]])
    - log10(len1 / len2)
}

go_graph$ict <- sapply(go_graph$term, get_ict)


#取ict值小于阈值的
nodes_cutoff <- mapply(function(e, d) {
    go_graph[go_graph$ont == e & go_graph$ict <= d, ]$term
}, c("b", "m", "c"), c(3.2, 3.6, 3.0))


#nodes_cutoff中存在某节点与其父节点们有过于close情况
#(ict比值小于1.2)
close_proximity <- function(terms, par = parents, g_graph = go_graph) {
    #事先预存好所有的节点
    all_ <- terms
    for (term1 in terms) {
        #父节点
        obj <- intersect(par[[term1]], all_)
        for (term2 in obj) {
            ict1 <- g_graph[g_graph[["term"]] == term1, ]$ict
            ict2 <- g_graph[g_graph[["term"]] == term2, ]$ict
            #判断
            if (ict2 != 0 & ict1 / ict2 <= 1.2) {
                #存在此情况即删除
                terms <- setdiff(terms, term1)
                break
            }
        }
    }
    #返回处理过之后的所有节点
    return(terms)
}


#得到各个ont去冗余之后的sub-root-nodes
meta_terms <- lapply(nodes_cutoff, close_proximity)

sub_root_nodes <- unlist(unname(meta_terms))


#以每个term的offspring作为meta_terms的内容，组成sub-graph
#每个sub-graph去掉来自其他sub-graph的信息：后代节点
remove_dup <- function(id, g_graph = go_graph, meta_term = sub_root_nodes) {
    content <- g_graph[g_graph[["term"]] == id, ]
    #该节点的所有后代
    offs <- unlist(content$offspring)
    #其他的sub-root-nodes
    other_sub <- intersect(unlist(offs), meta_term)
    #去掉本身
    other_sub <- setdiff(other_sub, id)
    #其他sub-root-nodes的后代
    loca <- charmatch(other_sub, g_graph[["term"]])
    other_sub_offs <- unlist(g_graph[loca, ]$offspring)
    #去掉
    unlist(setdiff(offs, other_sub_offs))
}





#meta_graph有3列，id，ont, terms
#分别代表:sub-root-id, ontology, sub-graph-terms
#存为data.frame格式，便于操作
meta_graph <- data.frame(id = sub_root_nodes,
                         ont = rep(c("b", "m", "c"),
                                   unlist(lapply(meta_terms, length))),
                         terms = I(lapply(sub_root_nodes, remove_dup)),
                         stringsAsFactors = F)


#创建"meta"作为一个cluster，将各个ont的sub-root-nodes集结到一起
tmp <- data.frame(id = rep("meta", 3),
                  ont = c("b", "m", "c"),
                  terms = I(meta_terms), stringsAsFactors = F)

#添加
meta_graph <- rbind(meta_graph, tmp)
