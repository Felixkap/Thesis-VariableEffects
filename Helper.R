library(tidyverse)

################################################################################
################################################################################
################################################################################
####################### True Relation between X_i and Y#########################
################################################################################
################################################################################
################################################################################

split.formula <- function(form, char_splits){
  preds <- strsplit(form, "(?<=[-+])", perl=TRUE)[[1]]
  print(preds)
  final_preds <- character()
  cnt <- 1 
  for (i in preds) {
    if(any(sapply(c('\\+', '\\-'), grepl, i))){
      final_preds[cnt] <- substring(i, 1, nchar(i)-1)
      print(final_preds[cnt])
      final_preds[cnt+1] <- substring(i, nchar(i), nchar(i))
      print(final_preds[cnt+1])
      cnt = cnt + 2
    } else {
      final_preds[cnt] <- i
      print(final_preds[cnt])
      cnt =+ 1
    }
  }
  return(final_preds)
}



inbetween_integers <- function(a, b) {
  u <- sort(c(a, b))
  res <- setdiff(ceiling(u[1]):floor(u[2]), c(a, b))
  if (!length(res)) {
    NULL
  } else {
    res
  }
}



fixed.formula.splitted <- split.formula(fixed.formula[2])

layout(matrix(c(1,2,
                3,4,
                5,5), nrow=3, byrow=T))
for (i in 1:5) {
  pred <- paste0("X.",i)
  idx <- grep(pred, fixed.formula.splitted)
  idx <- ifelse(is.null(inbetween_integers(idx[1], idx[length(idx)])), idx, c(idx[1], inbetween_integers(idx[1], idx[2]), idx[2]))
  
  form <- paste(fixed.formula.splitted[idx-1],fixed.formula.splitted[idx], sep='', collapse='')
  print(form)
  form <- gsub(pred, 'x', form)
  plot(function(x) eval(parse(text=form)), xlim = c(-3,3), ylab='Y', xlab=pred)
}



