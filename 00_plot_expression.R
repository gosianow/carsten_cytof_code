


plot_expression <- function(ggdf, color_groups, color_groupsb, outdir, prefix){
  # ggdf has to have the following columns: cluster, expr, group, data, day, marker
  

  # --------------------------------------------------------------------------
  # Plot each cluster as a separate page in the pdf file
  # --------------------------------------------------------------------------
  
  clusters <- levels(ggdf$cluster)
  
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
  
  
  
  # --------------------------------------------------------------------------
  # Plot each marker as a separate pdf file
  # --------------------------------------------------------------------------
  
  
  markers <- levels(ggdf$marker)
  days <- levels(ggdf$day)
  
  nr_cluster <- nlevels(ggdf$cluster)
  nr_in_one_row <- 5
  nrow <- ceiling(nr_cluster/nr_in_one_row)
  hh <- nrow * 2.5
  if(nr_cluster < nr_in_one_row){
    ww <- nr_cluster * 2 + 1
  }else{
    ww <- 2 * nr_in_one_row + 1
  }
  
  for(j in 1:nlevels(ggdf$marker)){
    # j = 1
    
    for(i in 1:nlevels(ggdf$day)){
      # j = 1; i = 1
      
      df <- ggdf[ggdf$marker == markers[j] & ggdf$day == days[i], , drop = FALSE]
      
      ggp <- ggplot(df) +
        geom_boxplot(aes(x = group, y = expr, color = group, fill = group), outlier.colour = NA, width = 0.9, position = position_dodge(width = 1)) +
        geom_point(aes(x = group, y = expr, color = group, shape = data), size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 1)) +
        ggtitle(markers[j]) +
        theme_bw() +
        ylab("Expression") +
        xlab("") +
        theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
          legend.title = element_blank(), legend.position = "right", legend.key = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0)) +
        scale_color_manual(values = color_groups) +
        scale_fill_manual(values = color_groupsb) +
        facet_wrap(~ cluster, scales = "free", ncol = nr_in_one_row)
      
      pdf(file.path(outdir, paste0(prefix, "expr_plot_marker_", gsub("-", "", markers[j]), "_", days[i] ,"2.pdf")), w = ww, h = hh)
      print(ggp)
      dev.off()
      
      
    }
    
    
    
  }
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  return(NULL)
  
}



