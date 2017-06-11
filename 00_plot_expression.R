


plot_expression <- function(ggdf, color_groups, color_groupsb, outdir, prefix){
  # ggdf has to have the following columns: cluster, expr, group, data, day, marker
  
  clusters <- levels(ggdf$cluster)
  
  # ------------------------------------
  ## plot each cluster as a separate page in the pdf file
  
  ggp <- list()
  
  for(i in 1:nlevels(ggdf$cluster)){
    # i = 1
    
    df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
    
    ggp[[i]] <- ggplot(df, aes(x = group, y = expr, color = group, shape = data)) +
    geom_boxplot(width = 0.9, position = position_dodge(width = 0.95), outlier.colour = NA) +
    geom_point(size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
      ggtitle(clusters[i]) +
      theme_bw() +
      ylab("Expression") +
      xlab("") +
      theme(axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "none") +
      scale_color_manual(values = color_groups) +
      facet_wrap(~ marker, scales = "free")
    
  }
  
  pdf(file.path(outdir, paste0(prefix, "expr_plot_both2.pdf")), w = 18, h = 12, onefile=TRUE)
  for(i in seq(length(ggp)))
    print(ggp[[i]])
  dev.off()
  
  
  return(NULL)
  
}



