setwd("D:/study/Yu/2020_stay_at_home/git_tcss/")

#处理源文件一，得到GO数据库节点之间的父子关系
#得到两个数据 parents children


#源数据一
#处理源数据一得到父子关系 ： parents
rr <- readLines("gene_ontology.obo.txt")

#取term所开始的行数
loca_s <- which(rr == "[Term]")

#id位于[term]下一行
id <- lapply(rr[loca_s + 1], function(e) sub("id: ", "", e))

#ontology位于[term]下三行
ont <- lapply(rr[loca_s + 3], function(e) substr(e, 12, 12))
#"b", "m", "c"

#所有的节点，包括id和ontology信息
total_node <- data.frame(id = unlist(id), ont = unlist(ont),
                         stringsAsFactors = F)


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
get_parents <- function(start, end, rr_ = rr) {
    lines <- rr_[start:end]
    unlist(lapply(lines, get_parent))
}

#得到每个term的parent terms，数据类型 list
#length : 33703, 是文件中所有的节点
#parents的索引是id，同时等于 total$id
parents <- mapply(get_parents, loca_s, loca_d)




#父子关系的进一步处理，子父关系 : children
#根据变量parents反向形成children
#我的想法：将parents展开，形成flip : term---child一一对应
rep_times <- unlist(lapply(parents, length))
rep_terms <- rep((as.character(id)), rep_times)
unfold_par <- unlist(parents)
flip <- data.frame(parent = unfold_par, child = rep_terms, stringsAsFactors = F)

#从flip中找到该term对应的所有child
#得到每个term的child terms，数据类型 list
#length : 20953 是所有节点中属于BP的节点
children <- lapply(total_node$id, function(e) {
    unique(flip[flip$parent == e, ]$child)
})

names(children) <- total_node$id
