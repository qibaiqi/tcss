#本文件旨在探究不同的cutoff下的模型预测PPI的效果
#主要通过画ROC曲线
#横轴是F1值，纵轴是曲线下面积

#源文件三
#已经确定为PPI的数据
#positive_set
#每一行是一对PPI，分别取第一个蛋白和第二个蛋白
r <- readLines("positives.human.iea_yes.p")
first_pro <- lapply(r, function(e) strsplit(e, ",")[[1]][1])
second_pro <- lapply(r, function(e) strsplit(e, ",")[[1]][2])
p_set <- data.frame("first_pro" = unlist(first_pro),
                    "second_pro" = unlist(second_pro),
                    stringsAsFactors = FALSE)

#源文件四
#negative_set
r <- readLines("negatives.human.iea_yes.p")
first_pro <- lapply(r, function(e) strsplit(e, ",")[[1]][1])
second_pro <- lapply(r, function(e) strsplit(e, ",")[[1]][2])
n_set <- data.frame("first_pro" = unlist(first_pro),
                    "second_pro" = unlist(second_pro),
                    stringsAsFactors = FALSE)

#add lable
label_pro <- c(rep("Yes", times = 1435), rep("No", times = 1435))#加标签
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
    #最佳阈值
    thres <- coords(roc1, "best", "threshold", transpose = FALSE)[1]
    thres <- as.numeric(thres)
    #根据阈值统计预测结果
    predict_label <- lapply(predict_value, function(e)
        if (is.nan(e)) "No" else (if (e > thres) "Yes" else "No"))
    predict_label <- unlist(predict_label)
    #得到混淆矩阵
    c_matrix <- table(predict_label, label_pro)
    #计算F1值
    2 * c_matrix[2, 2] / (2 * c_matrix[2, 2] + c_matrix[1, 2] + c_matrix[2, 1])
}


#运行
s <- Sys.time()
#用函数计算这些蛋白对
predict_value <- mapply(protcss, test_set$first_pro, test_set$second_pro)
#用roc计算不同阈值下的敏感度，精确度
#predict_value有空值，NULL，所以此处需要先替换一下
predict_value <- lapply(predict_value, function(e) if (is.null(e)) NaN else e)
predict_value <- unlist(predict_value)
roc1 <- roc(label_pro, predict_value, print.thres = TRUE, legacy.axes = TRUE)
#同时得到曲线下面积
auc <- as.numeric(auc(roc1))
#计算F1值
f1_score <- f1_score_calculate(roc1, predict_value)
e <- Sys.time()
print(e - s)
