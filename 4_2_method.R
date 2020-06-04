#蛋白对应多个term，term分为不同的cluster，产生多个value
#以value为矩阵或向量的形式传入，根据不同计算方式得到最终值



operate_method <- function(value, method = "max") {
    
    value <- unlist(value)
    if (method == "max") {
        #result <- max(unlist(value))
        result <- max(value, na.rm = TRUE)
    }
    if (method == "avg") {
        result <- mean(value, na.rm = TRUE)
    }
    if (method == "bma") {
        if (is.matrix(value)) {
            result <- sum( apply(value, 1, max, na.rm = TRUE),
                           apply(value, 2, max, na.rm = TRUE)
            ) / sum(dim(value))
        } else {
            result <- mean(value, na.rm = TRUE)
        }
    
    }
    return(round(result, 8))
}