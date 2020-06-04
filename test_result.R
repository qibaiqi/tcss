#此文件只为测试结果是否一致
#源文件五
#已经确定为PPI的数据
gg <- readLines("species human_MINT_PPI")
gg <- strsplit(gg, split = "\t")

#每一行是一对PPI，分别取第一个蛋白和第二个蛋白
#切的是不同列
first_pro <- unlist(lapply(gg, function(e) strsplit(e[1], ":")[[1]][2]))
second_pro <- unlist(lapply(gg, function(e) strsplit(e[2], ":")[[1]][2]))

#因为有些protein在goa_human.gaf中并没有被注释到,而去除
#老师的方法 用length判断是否存在
#组成PPI组合
ppi_pair <- data.frame(first_pro, second_pro, stringsAsFactors = F)
len1 <- sapply(pro_annotations[ppi_pair$first_pro], length)
len2 <- sapply(pro_annotations[ppi_pair$second_pro], length)
ppi_pair <- ppi_pair[len1 & len2, ]
#pair_set去重 保留了protein自己与自己本身,20637行
ppi_pair <- ppi_pair[!duplicated(ppi_pair), ]

#共10000对PPI，由正向5000和负向5000组成
test_sets <- rbind(ppi_pair[1:5000, ],
                   data.frame(first_pro = sample(first_pro, 5000, replace = T),
                       second_pro = sample(second_pro, 5000, replace = T)))


run_with_python <- paste0(test_sets$first_pro, ",", test_sets$second_pro)
write.table(run_with_python, "test_file.txt", quote = F, row.names = FALSE,
            col.names = FALSE)

##本地python2.6+运行
#python2 tcss.py --gene gene_association.goa_human -c P:3.2 -i test_file.txt -o result_file.txt

##本地目录手动复制
#读取作者的结果
tt <- readLines("result_file.txt")
tt <- strsplit(tt, split = "\t")


get_author_value <- function(line) {
    if (length(line) == 0) {
        return(NULL)
    }else if (grepl("Biological", line) |
              grepl("Molecular", line) |
              grepl("Cellular", line)) {
        value <- as.numeric(sub("\\D+", "\\1", line))
        return(value)
    }
}

author_result <- unlist(lapply(tt, get_author_value))
author_result <- round(author_result, 8)

#我的结果
#cutoff须保持一致
my_result <- mapply(protcss, test_sets$first_pro, test_sets$second_pro, "b")
#my_result <- mapply(protcss, test_sets$first_pro, test_sets$second_pro, "m")
#my_result <- mapply(protcss, test_sets$first_pro, test_sets$second_pro, "c")
my_result <- lapply(my_result, function(e) if (is.null(e)) NA else e)
my_result <- unlist(my_result)
my_result <- round(my_result, 8)
my_result <- unname(my_result)

combine <- mapply(function(e, f) all.equal(e, f), author_result, my_result)
table(combine)
identical(author_result, my_result)
