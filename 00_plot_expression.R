


plot_expression <- function(ggdf, ggds, color_groups, outdir, prefix, prefix2){
  
  
  clusters <- levels(ggdf$cluster)
  
  # ------------------------------------
  ## plot each cluster as a separate page in the pdf file
  
  ggp <- list()
  
  for(i in 1:nlevels(ggdf$cluster)){
    # i = 1
    
    df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
    ds <- ggds[ggds$cluster == clusters[i], , drop = FALSE]
    
    ggp[[i]] <- ggplot(df, aes(x = group, y = expr, color = group, shape = data)) +
      geom_point(size=2, position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.7)) +
      geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean, ymax=mean), color='black', width=0.4, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
      geom_errorbar(data=ds, aes(x=group, y=mean, ymin=mean-sd, ymax=mean+sd), color='black', width=0.25, position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.7)) +
      facet_wrap(~ marker, scales = "free") +
      ggtitle(clusters[i]) +
      theme_bw() +
      ylab("Expression") +
      xlab("") +
      theme(axis.text.x = element_text(size=10, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "right") +
      scale_color_manual(values = color_groups)
    
  }
  
  pdf(file.path(outdir, paste0(prefix, "expr_", prefix2, ".pdf")), w = 18, h = 12, onefile=TRUE)
  for(i in seq(length(ggp)))
    print(ggp[[i]])
  dev.off()
  
  
  return(NULL)
  
}



