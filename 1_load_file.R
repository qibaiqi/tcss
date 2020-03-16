setwd("D:/Pictures/important files/2020_stay_at_home/git_tcss/")
#源数据一
rr <- readLines("gene_ontology.obo.txt")
loca_s <- which(rr == "[Term]")#取term所开始的行数

#提取term的id
get_id <- function(line) {
  if (grepl("id", line)) {
    id <- sub("id: ", "", line)
    return(id)
  }
}

#提取term的ontology
get_ont <- function(line) {
  if (grepl("namespace", line)) {
    ont <- substr(line, 12, 12)
    return(ont)
  }
}

id <- lapply(rr[loca_s + 1], get_id)#id位于[term]下一行
ont <- lapply(rr[loca_s + 3], get_ont)#ontology位于[term]下三行

node <- data.frame(id = unlist(id), ont = unlist(ont), stringsAsFactors = F)
node_bp <- node[node$ont == "b", ]$id#只选取ontology为BP的node

loca_d <- c(rep(loca_s, 1), length(rr))[-1]#每个term的信息结束的行数

#提取term的parent列表
#get_parent处理一行(one line)的信息
get_parent <- function(line) {
  if (grepl("is_a", line) | grepl("relationship", line)) {
    parent <- sub(".*(GO:(\\d+)).*", "\\1", line)
    return(parent)
  }else {
      return(NULL)
    }
}

#get_parents处理[loca_s, loca_d]一个段落(several lines)的信息
get_parents <- function(start, end) {
  lines <- rr[start:end]
  parents <- unlist(lapply(lines, get_parent))
  return(parents)
}

parents <- mapply(get_parents, loca_s, loca_d)
names(parents) <- id


#parents累加形成ancestors
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
  return(items)
}

ancestors <- lapply(node_bp, depth_search, parents)
names(ancestors) <- node_bp


#根据变量parents反向形成children
#我的想法：将parents展开，形成flip : term---child一一对应
rep_times <- unname(unlist(lapply(parents, length)))
rep_terms <- rep((as.character(id)), rep_times)
unfold_par <- unlist(unname(parents))
flip <- data.frame(parent = unfold_par, child = rep_terms, stringsAsFactors = F)

#从信息池(flip)中找到该term对应的所有child
get_children <- function(term, pool) {
  children <- pool[pool$parent == term, ]
  children <- unique(children$child)
  return(children)
}

children <- lapply(node_bp, get_children, flip)
names(children) <- node_bp

#children的所有children...迭代得到 offspring
offspring <- lapply(node_bp, depth_search, children)
names(offspring) <- node_bp

#源数据二
#term与gene的关联信息文件
ff <- readLines("gene_association.sgd")
ff <- ff[-c(1:28)]## 去掉前几行介绍信息
ff <- strsplit(ff, split = "\t")
ff <- t(as.data.frame(ff))## 变成了matrix
#文件第二列是gene id, 第五列是对应的 term id
gene_term <- data.frame(gene = ff[, 2], term = ff[, 5], stringsAsFactors = F)

#提取term的注释信息：gene列表
get_anno <- function(term, pool) {
  genes <- pool[pool$term == term, ]
  genes <- unique(genes$gene)
  return(genes)
}

#初始的annotations，每个term对应的gene列表
init_annotations <- lapply(node_bp, get_anno, gene_term)
names(init_annotations) <- node_bp

#annotations信息的合并
#将该term的所有后代(offspring)的注释信息全都加到自己身上
anno_add <- function(term) {
  off_list <- offspring[[term]]
  anno_list <- lapply(off_list, function(e) init_annotations[[e]])
  anno_list <- c(unlist(anno_list), init_annotations[[term]])
  anno_list <- unique(anno_list)
  return(anno_list)
}

#最终的annotations
final_annotations <- lapply(node_bp, anno_add)
names(final_annotations) <- node_bp

#从parents,children,ancestor,offspring,final_annotations
#中去掉注释为空的的term, delete

#final_annotations
dele_loca <- which(lapply(final_annotations, length) == 0)
node_remove <- node_bp[dele_loca]

final_annotations <- final_annotations[-dele_loca]


#以变量为对象，这样才能持续的改变变量
#delete_inside是从list的各个字符串中删除这些节点
delete_inside <- function(term_set, node_remove) {
  cut_terms <- intersect(term_set, node_remove)
  if (length(cut_terms) != 0) {
    term_set <- setdiff(term_set, cut_terms)
    return(term_set)
  }else{
    return(term_set)
  }
}

#children
children <- lapply(children, delete_inside, node_remove)
children[node_remove] <- NULL

#ancestors
ancestors <- lapply(ancestors, delete_inside, node_remove)
ancestors[node_remove] <- NULL

#offspring
offspring <- lapply(offspring, delete_inside, node_remove)
offspring[node_remove] <- NULL
