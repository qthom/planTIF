
findMatchedControl_v2 <- function(treat, ctrl, log=TRUE, pseudocount=1) {
  if (isTRUE(log)) {
    treat <- log2(treat + pseudocount)
    ctrl <- log2(ctrl + pseudocount)
  }
  breaks <- seq(min(ctrl), max(ctrl), length.out=100)
  treat_cut <- cut(treat, breaks=breaks)
  ctrl_cut <- cut(ctrl, breaks=breaks)
  df_treat <- as.data.frame(table(treat_cut))
  names(df_treat) <- c("bin", "freq")
  df_treat$bin <- as.numeric(df_treat$bin)
  df_treat <- df_treat[df_treat$freq > 0, ]
  df_ctrl <- data.frame("value"=ctrl, "index"=1:length(ctrl), "bin"=as.numeric(ctrl_cut))
  df_ctrl$bin <- as.numeric(df_ctrl$bin)
  df_ctrl <- df_ctrl[!is.na(df_ctrl$bin), ]
  results <- vector("list", nrow(df_treat))
  for (i in 1:nrow(df_treat)) {
    freq <- df_treat$freq[i]
    bin <- df_treat$bin[i]
    indexes <- df_ctrl$index[df_ctrl$bin==bin]
    freq <- min(freq, length(indexes))
    results[[i]] <- sample(indexes, size=freq)
  }
  out <- sort(do.call(c, results))
  return(out)
}
