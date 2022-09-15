#Bootstrap confidence interval fuction
conf_int <- function(vector, level = 0.95, boot.iter = 1000){
  bootstrap = sapply(1:boot.iter, function(x){
    mean(sample(vector, replace = T))
  })
  lower <- quantile(bootstrap, c((1-level) / 2))
  upper <- quantile(bootstrap, c(level + (1-level) / 2))
  
  c(upper = upper, lower = lower)
}

#p-value sigificance symbols
sig_symbols <- function(p_vals, non_sig_symbol = "N.S"){
  
  #takes a numeric vector of p-values as input and outputs a vector significance
  #symbols
  
  len <- length(p_vals)
  symbol <- rep(non_sig_symbol, len)
  
  for (i in 1:len){
    
    if (p_vals[i] < 0.001){
      symbol[i] <- "***"
    } else if (p_vals[i] < 0.01){
      symbol[i] <- "**"
    } else if (p_vals[i] < 0.05){
      symbol[i] <- "*"
    }
  } 
  return(symbol)
}

#cohen's D
cohens_d = function(x,y){ 
  
  #calculates effect size of vectors a and b, using cohen's d by http://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d
  lx=length(x)
  ly=length(y)
  (mean(x)-mean(y)) / ( ( ( (lx-1)*var(x) + (ly-1)*var(y) ) / ( lx+ly-2 ) )^0.5 )
}

#Subsetting function
subsetCommon <- function(matrix1, colname1 = NULL, matrix2, colname2 = NULL){
  #A method to subset two different matrix-like objects based on column values
  #if no columns are given, rownames are used by default. 
  #Outputs is a list of ordered, subsetted matrix-like objects
  
  if (is.null(colname1) && is.null(colname2)){
    matrix1_subset <- matrix1[rownames(matrix1) %in% rownames(matrix2), ]
    matrix1_subset <- matrix1_subset[order(rownames(matrix1_subset)), ]
    
    matrix2_subset <- matrix2[rownames(matrix2) %in% rownames(matrix1), ]
    matrix2_subset <- matrix2_subset[order(rownames(matrix2_subset)), ]
  } else if (is.null(colname1) && !is.null(colname2)) {
    
    matrix1_subset <- matrix1[rownames(matrix1) %in% matrix2[,colname2], ]
    matrix1_subset <- matrix1_subset[order(rownames(matrix1_subset)), ]
    
    matrix2_subset <- matrix2[matrix2[,colname2] %in% rownames(matrix1), ]
    matrix2_subset <- matrix2_subset[order(matrix2_subset[,colname2]), ]
    
  } else if (!is.null(colname1) && is.null(colname2)) {
    
    matrix1_subset <- matrix1[matrix1[,colname1] %in% rownames(matrix2), ]
    matrix1_subset <- matrix1_subset[order(matrix1_subset[,colname1]), ]
    
    matrix2_subset <- matrix2[rownames(matrix2) %in% matrix1[,colname1], ]
    matrix2_subset <- matrix2_subset[order(rownames(matrix2_subset)), ]
    
  } else if (!is.null(colname1) && !is.null(colname2)) {
    
    matrix1_subset <- matrix1[matrix1[,colname1] %in% matrix2[,colname2], ]
    matrix1_subset <- matrix1_subset[order(matrix1_subset[,colname1]), ]
    
    matrix2_subset <- matrix2[rownames(matrix2) %in% matrix1[,colname1], ]
    matrix2_subset <- matrix2_subset[order(matrix2_subset[,colname2]), ]
    
  }
  
  a_list <- list(matrix1 = matrix1_subset, matrix2 = matrix2_subset)
  a_list
}

