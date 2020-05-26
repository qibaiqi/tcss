#处理源数据二,得到注释信息:init_annotations, pro_term
#源数据二
#term与gene的关联信息文件



#此为酵母蛋白注释信息
# ff <- readLines("gene_association.sgd")
# ff <- ff[-c(1:28)]## 去掉前几行介绍信息
# ff <- strsplit(ff, split = "\t")
# ff <- t(as.data.frame(ff))## 变成了matrix
# #文件第二列是gene id, 第五列是对应的 term id
# pro_term <- data.frame(pro = ff[, 2], term = ff[, 5], stringsAsFactors = F)



#此为人类蛋白注释信息
dd <- readLines("gene_association.goa_human")
dd <- dd[!grepl("!", dd)]## 注释行以"!"开头
dd <- strsplit(dd, split = "\t")
# 变成了matrix,前1447行是17列，后面是16列
ff <- t(as.data.frame(dd[1:1447]))
ff <- ff[, 1:16]
ff <- rbind(ff, t(as.data.frame(dd[1448:length(dd)])))
#文件第二列:gene id, 第五列:term id, 第三列:evidence code
pro_term <- data.frame(pro = ff[, 2], term = ff[, 5], stringsAsFactors = F)
#pro_term <- data.frame(pro = ff[, 2], term = ff[, 5], code = ff[, 7],
#                       stringsAsFactors = F)



#初始的annotations，每个term对应的gene列表
init_annotations <- lapply(total_node$id, function(e) {
    unique(pro_term[pro_term$term == e, ]$pro)
})
names(init_annotations) <- total_node$id
