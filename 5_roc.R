#本文件旨在探究不同的cutoff下的模型预测PPI的效果，主要通过化ROC曲线
#横轴是F1值，纵轴是曲线下面积

#源文件三
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
len1 <- sapply(gene_annotations[ppi_pair$first_pro], length)
len2 <- sapply(gene_annotations[ppi_pair$second_pro], length)
ppi_pair <- ppi_pair[len1 & len2, ]
#pair_set去重 保留了protein自己与自己本身
ppi_pair <- ppi_pair[!duplicated(ppi_pair), ]


#positive_set
p_set <- ppi_pair[1:10000, ]#选取1000行

#negative_set
#随机选取1000对
n_set <- data.frame(first_pro = sample(ppi_pair$first_pro, 10000, replace = T),
                    second_pro = sample(ppi_pair$second_pro, 1000, replace = T))

#add lable
label_pro <- c(rep("Yes", times = 10000), rep("No", times = 10000))#加标签
test_set <- rbind(p_set, n_set)


#install.packages("pROC")
#library("pROC")
#data(aSAH)
#aSAH
#roc1<-roc(aSAH$outcome, aSAH$s100b, plot = TRUE,
#print.thres = TRUE,legacy.axes=TRUE)
#                label_pro
# predict_label   No  Yes
#            No  2501 4646
#            Yes 7499 5354


#根据roc的结果和预测结果计算F1值
f1_score_calculate <- function(roc1, predict_value) {
    #根据约登指数选出最佳阈值
    yuoden <- mapply(function(x, y) x + y - 1,
                       roc1$sensitivities, roc1$specificities)
    #得到最佳阈值的位置
    pos <- charmatch(max(yuoden), yuoden)
    #得到最佳阈值
    thres <- roc1$thresholds[pos]
    #根据阈值统计预测结果
    predict_label <- lapply(predict_value,
                            function(e) if (e < thres) "Yes" else "No")
    predict_label <- unlist(predict_label)
    #得到混淆矩阵
    c_matrix <- table(predict_label, label_pro)
    #计算F1值
    2 * c_matrix[2, 2] / (2 * c_matrix[2, 2] + c_matrix[1, 2] + c_matrix[2, 1])
}


#运行
s <- Sys.time()
#用函数计算这一万个蛋白对
predict_value <- mapply(protcss, test_set$first_pro, test_set$second_pro)
#用roc计算不同阈值下的敏感度，精确度
roc1 <- roc(label_pro, predict_value, plot = FALSE,
                print.thres = TRUE, legacy.axes = TRUE)
#同时得到曲线下面积
auc <- roc1$auc[1]
#计算F1值
f1_score <- f1_score_calculate(roc1, predict_value)
e <- Sys.time()
print(e - s)
