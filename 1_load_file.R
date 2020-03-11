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
get_anno <- function(term, pool) {
  genes <- pool[pool$term == term, ]
  genes <- unique(genes$gene)
  return(genes)
}

annotations <- lapply(node_bp, get_anno, gene_term)
names(annotations) <- node_bp


#根据变量parents反向形成children
rep_times <- unname(unlist(lapply(parents, length)))
rep_terms <- rep((as.character(id)), rep_times)
unfold_par <- unlist(unname(parents))
flip <- data.frame(parent = unfold_par, child = rep_terms, stringsAsFactors = F)

get_children <- function(term, pool) {
  children <- pool[pool$parent == term, ]
  children <- unique(children$child)
  return(children)
}

children <- lapply(id, get_children, flip)
names(children) <- id
