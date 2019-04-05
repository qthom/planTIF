processCoords <- function(gr, signal, expand) {
  expanded_gr <- trim(resize(gr, expand, fix="center"))
  expanded_gr <- expanded_gr[width(expanded_gr) == expand]
  mat <- metageneMatrix(signal, expanded_gr, scale=FALSE, skip.outliers=FALSE)
  avg <- apply(mat, 2, mean)
  sem <- apply(mat, 2, sd) / sqrt(nrow(mat))
  lower <- avg - 2 * sem
  upper <- avg + 2 * sem
  out <- list("avg" = avg, "lower" = lower, "upper" = upper)
}

parseResults <- function(x, y, expand = expand) {
  df <- as.data.frame(do.call(cbind, x))
  df$group <- y
  df$pos <- seq(-expand/2, length.out = expand)
  df
}

metagenePipeline <- function(signal, intervals, skip=NULL, out_dir=".", expand) {
  output <- list()
  for (i in seq_along(signal)) {
    sample_name <- names(signal)[[i]]
    if (!is.null(skip)) {
      sample_name <- gsub(skip, "", sample_name)
    }
    cat(sample_name, "\n"); flush.console()
    data <- signal[[i]]
    results <- lapply(intervals, processCoords, signal=data, expand=expand)
    cat("parseResults working...\n"); flush.console()
    parsed <- mapply(results, names(results), FUN=parseResults, SIMPLIFY=F, expand = expand)
    cat("Making dataframe...\n")
    df <- do.call(rbind, parsed)
    df <- subset(df, lower >= 0)
    suppressWarnings(dir.create(out_dir))
    out_name <- file.path(out_dir, sample_name)
    cat("Plotting...\n")
    cbPalette <- c("#0072B2","darkorange2")
    p <- ggplot(df, aes(x = pos, y = avg, color = group)) + theme_classic() + geom_line(size=2) + geom_vline(xintercept = 0) + geom_ribbon(aes(x = pos, ymin = lower, ymax = upper, fill = group),linetype = "blank", alpha = 0.45) + ylim(0, NA)  + ggtitle(sample_name) + xlab("Relative position from TSS") + ylab(yaxlab) + scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
  
    ggsave(filename = paste0(out_name, ".png"), plot = p)
    ggsave(filename = paste0(out_name, ".pdf"), plot = p)
    cat("Figure was saved;\n\n")
    output[[i]] <- sample_name
  }
  output
}



theme_classic <- function(base_size = 11, base_family = "",
                          base_line_size = base_size / 22,
                          base_rect_size = base_size / 22) {
  theme_bw(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    theme(
      # no background and no grid
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # show axes
      axis.line      = element_line(colour = "black", size = rel(1)),

      # match legend key to panel.background
      legend.position       = "bottom",

      # simple, black and white strips
      strip.background = element_rect(fill = "white", colour = "black", size = rel(2)),
      # NB: size is 1 but clipped, it looks like the 0.5 of the axes

      complete = TRUE
    )
}
