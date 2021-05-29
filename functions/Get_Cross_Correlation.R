#Function to get a plot of the Cross Correlation 

#Input
#1. Dataset
#2. Title of Plot



get_crosscorrelations <- function(input_data,nam){


#Load the Input Data
x <- input_data
cols <- ncol(x)
  
    
  
###Load the Packages
library(corrplot)


###Computing the Correlations
M<-cor(x)


# mat : is a matrix of data. Source:- http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# matrix of the p-value of the correlation
p.mat <- cor.mtest(x)


#Small Datasets
if(cols < 99){
  corrplot(M, type="upper", method = "shade",
           p.mat = p.mat, sig.level = 0.01, insig = "blank", 
          title = paste0("Cross Correlation \n ", nam), tl.pos = "n",
           mar = c(1,1,4,1))
}




#Large Datasets
if(cols > 99){
corrplot(M, type="upper", order="hclust", method = "shade",
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         tl.pos = "n", title = paste0("Cross Correlation  \n ", nam), 
         mar = c(1,1,4,1))
  }

}