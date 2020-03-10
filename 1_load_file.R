setwd("D:/Pictures/important files/2020_stay_at_home/git_tcss/")
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

id <- lapply(rr[loca_s + 1], get_id)
ont <- lapply(rr[loca_s + 3], get_ont)

node <- data.frame(id = unlist(id), ont = unlist(ont), stringsAsFactors = F)
node_bp <- node[node$ont == "b", ]$id#只选取BP

loca_d <- c(rep(loca_s, 1), length(rr))[-1]#每个term的信息结束的行数

#提取term的parent列表
get_parent <- function(line) {
  if (grepl("is_a", line) | grepl("relationship", line)) {
    parent <- sub(".*(GO:(\\d+)).*", "\\1", line)
    return(parent)
  }else {
      return(NULL)
    }
}

get_parents <- function(start, end) {
  lines <- rr[start:end]
  parents <- unlist(lapply(lines, get_parent))
  return(parents)
}

parents <- mapply(get_parents, loca_s, loca_d)
names(parents) <- id

#term与gene的关联信息文件
ff <- readLines("gene_association.sgd")
ff <- ff[-c(1:28)]## 去掉前几行介绍信息
ff <- strsplit(ff, split = "\t")
ff <- t(as.data.frame(ff))## 变成了matrix
gene_term <- data.frame(gene = ff[, 2], term = ff[, 5], stringsAsFactors = F)

#提取term的注释信息：gene列表
get_anno <- function(term) {
  genes <- gene_term[gene_term$term == term, ]
  genes <- unique(genes$gene)
  return(genes)
}

annotations <- lapply(node_bp, get_anno)
names(annotations) <- node_bp


save(list = ls(), file = "1_load_file.RData")

#注释信息为空的行数
term_remove_loca <- which(lapply(annotations, length) == 0)


###########
annotations <- annotations[-term_remove_loca]

term_remove <- node_bp[term_remove_loca]

node_bp <- node_bp[-term_remove_loca]

parents <- lapply(term_remove, function(e) parents[[e]] <- NULL)

#flip
rep_times <- unname(unlist(lapply(parents, length)))
rep_terms <- rep()


###########
rep_times <- unname(unlist(lapply(ff, length)))
rep_terms <- rep(term_list, rep_times)
flip_child <- data.frame(term = children, child = rep_terms)

get_ont_children <- function(item) {
  child <- subset(flip_child, flip_child$term == item, select = c("child"))
  child <- as.character(unique(child)$child)
}
ont_children <- lapply(term_list, get_ont_children)
names(ont_children) <- term_list
ont_children <- ont_children[lapply(ont_children, length) > 0]








#[Term]
#id: GO:0000018
#name: regulation of DNA recombination
#namespace: biological_process
#def: "Any process that modulates ...
#subset: gosubset_prok
#is_a: GO:0051052 ! regulation of DNA metabolic process
#relationship: regulates GO:0006310 ! DNA recombination
