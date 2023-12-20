


#' @description softmax function
#' @param x a numeric vector
#' @return a normalized numeric vector
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @details
#' 尽管softmax是一个非线性函数，但softmax回归的输出仍然由输入特征的仿射变换决定。因此，softmax回归是一个线性模型（linear model）。
#' @examples
#' a <- softmax(c(1,2,3,4))
#' sum(a)
softmax <- function(x){
  return(exp(x)/sum(exp(x)))
}
