setwd("D:/study/Yu/2020_stay_at_home/git_tcss/")

#本文件旨在处理两份源文件，第一份得到GO数据库所有term之间的父子关系，
#并整理为children, parents, ancestors, offspring四种，第二份源文件，
#整理得到term和protein之间的关系，并形成注释，此处得到term的注释

#源数据一
rr <- readLines("gene_ontology.obo.txt")

#取term所开始的行数
loca_s <- which(rr == "[Term]")

#id位于[term]下一行
id <- lapply(rr[loca_s + 1], function(e) sub("id: ", "", e))

#ontology位于[term]下三行
ont <- lapply(rr[loca_s + 3], function(e) substr(e, 12, 12))

#所有的节点，包括id和ontology信息
node <- data.frame(id = unlist(id), ont = unlist(ont), stringsAsFactors = F)
##只选取ontology为BP的node
node_bp <- node[node$ont == "b", ]$id

#每个term的信息结束的行数
loca_d <- c(loca_s, length(rr))[-1]

#提取term的parent列表
#get_parent处理一行(one line)的信息
get_parent <- function(line) {
    if (grepl("is_a:", line) | grepl("relationship:", line)) {
        sub(".*(GO:(\\d+)).*", "\\1", line)
        }
    }

#get_parents处理[loca_s, loca_d]一个段落(several lines)的信息
get_parents <- function(start, end) {
    lines <- rr[start:end]
    unlist(lapply(lines, get_parent))
}

#得到每个term的parent terms，数据类型 list
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
get_children <- function(term, pool = flip) {
    part_pool <- pool[pool$parent == term, ]
    unique(part_pool$child)
}

children <- lapply(node_bp, get_children)
names(children) <- node_bp

#children的所有children...迭代得到 offspring
offspring <- lapply(node_bp, depth_search, children)
names(offspring) <- node_bp

#源数据二
#term与gene的关联信息文件
# ff <- readLines("gene_association.sgd")
# ff <- ff[-c(1:28)]## 去掉前几行介绍信息
# ff <- strsplit(ff, split = "\t")
# ff <- t(as.data.frame(ff))## 变成了matrix
# #文件第二列是gene id, 第五列是对应的 term id
# gene_term <- data.frame(gene = ff[, 2], term = ff[, 5], stringsAsFactors = F)

#此为人类蛋白注释信息，为方便后面的ROC
dd <- readLines("gene_association.goa_human")
dd <- dd[-c(1:25)]## 去掉前几行介绍信息
dd <- strsplit(dd, split = "\t")
## 变成了matrix,前1447行是17列，后面是16列
ff <- t(as.data.frame(dd[1:1447]))
ff <- ff[, 1:16]
ff <- rbind(ff, t(as.data.frame(dd[1448:192807])))
#文件第二列是gene id, 第五列是对应的 term id
pro_term <- data.frame(pro = ff[, 2], term = ff[, 5], stringsAsFactors = F)

#提取term的注释信息：gene列表
get_anno <- function(term, pro_terms = pro_term) {
    part_pro_terms <- pro_terms[pro_terms$term == term, ]
    unique(part_pro_terms$pro)
}

#初始的annotations，每个term对应的gene列表
init_annotations <- lapply(node_bp, get_anno)
names(init_annotations) <- node_bp

#annotations信息的合并
#父节点的注释信息等于其本身的注释信息及其所有后代节点的注释信息的总和
anno_add <- function(term, offs = offspring, annotations = init_annotations) {
    off_list <- offs[[term]]
    anno_list <- lapply(off_list, function(e) annotations[[e]])
    anno_list <- c(unlist(anno_list), annotations[[term]])
    anno_list <- unique(anno_list)
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
delete_inside <- function(term_set, remove = node_remove) {
    term_set <- setdiff(term_set, remove)
}

#children
children <- lapply(children, delete_inside)
children[node_remove] <- NULL

#ancestors
ancestors <- lapply(ancestors, delete_inside)
ancestors[node_remove] <- NULL

#offspring
offspring <- lapply(offspring, delete_inside)
offspring[node_remove] <- NULL

source("test_load_file.R")
source("2_clustering.R")
