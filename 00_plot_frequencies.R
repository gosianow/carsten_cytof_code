


plot_frequencies <- function(ggdf, color_groups, color_groupsb, colors_clusters = NULL, outdir, prefix){
  # ggdf has to have the following columns: cluster, prop, group, data, day
  
  if(is.null(ggdf)){
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot_both2.pdf")))
    plot(1, type="n", axes=F, xlab="", ylab="")
    dev.off()
    
  }else{
    # nr_cluster <- nlevels(ggdf$cluster)
    # 
    # ggp <- ggplot(ggdf, aes(x = cluster, y = prop, color = group, shape = data, fill = group)) +
    #   geom_boxplot(width = 0.95, position = position_dodge(width = 0.95), outlier.colour = NA) +
    #   geom_point(size=1.5, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
    #   theme_bw() +
    #   ylab("Frequency") +
    #   xlab("") +
    #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), 
    #     axis.title.y = element_text(size = 12), 
    #     panel.grid.major = element_blank(), 
    #     panel.grid.minor = element_blank(), 
    #     panel.border = element_blank(), 
    #     axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
    #     axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
    #     legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
    #   guides(color = guide_legend(ncol = 1)) +
    #   scale_color_manual(values = color_groups) +
    #   scale_fill_manual(values = color_groupsb) +
    #   facet_wrap(~ day)
    # 
    # pdf(file.path(outdir, paste0(prefix, "frequencies_plot_both0.pdf")), w = nr_cluster + 3, h = 4)
    # print(ggp)
    # dev.off()
    
    
    # ------------------------------------
    ### facet per cluster
    
    nr_cluster <- nlevels(ggdf$cluster)
    nr_in_one_row <- 5
    nrow <- ceiling(nr_cluster/nr_in_one_row)
    hh <- nrow * 2.5
    if(nr_cluster < nr_in_one_row){
      ww <- nr_cluster * 3 + 1
    }else{
      ww <- 3 * nr_in_one_row + 1
    }
    
    clusters <- levels(ggdf$cluster)
    
    ggp <- ggplot(ggdf, aes(x = group, y = prop, color = group, shape = data, fill = group)) +
      geom_boxplot(width = 0.9, position = position_dodge(width = 0.95), outlier.colour = NA) +
      geom_point(size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
      theme_bw() +
      ylab("Frequency") +
      xlab("") +
      theme(axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.position = "right", legend.key = element_blank(), legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0)) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = color_groups) +
      scale_fill_manual(values = color_groupsb) +
      facet_wrap(~ cluster, scales = "free", ncol = 4)
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot_both2.pdf")), w = ww, h = hh)
    print(ggp)
    dev.off()
    
    
    # ------------------------------------
    # plot all clusters in one pdf; colors per response; separate pdf for base and tx; one boxplot + points
    
    
    # days <- levels(ggdf$day)
    # nr_cluster <- nlevels(ggdf$cluster)
    # hh <- 4
    # ww = nr_cluster + 1
    
    # for(i in 1:nlevels(ggdf$day)){
    #   # i = 1
      
    #   df <- ggdf[ggdf$day == days[i], , drop = FALSE]
      
    #   ggp <- ggplot(df) +
    #     geom_boxplot(aes(x = cluster, y = prop, color = group, fill = group), width = 0.95, position = position_dodge(width = 0.95), outlier.colour = NA) +
    #     geom_point(aes(x = cluster, y = prop, color = group, shape = data), size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
    #     theme_bw() +
    #     ylab("Frequency") +
    #     xlab("") +
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12), 
    #       axis.title.y = element_text(size=12), 
    #       panel.grid.major = element_blank(), 
    #       panel.grid.minor = element_blank(), 
    #       panel.border = element_blank(), 
    #       axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
    #       axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
    #       legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
    #     scale_color_manual(values = color_groups) +
    #     scale_fill_manual(values = color_groupsb)
      
    #   pdf(file.path(outdir, paste0(prefix, "frequencies_plot_", days[i] ,"0.pdf")), w = ww, h = hh)
    #   print(ggp)
    #   dev.off()
      
    # }
    
    
    # ------------------------------------
    # plot all clusters in one pdf; colors per response; separate pdf for base and tx; one boxplot + points + facet per cluster
    
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
    
    for(i in 1:nlevels(ggdf$day)){
      # i = 1
      
      df <- ggdf[ggdf$day == days[i], , drop = FALSE]
      
      ggp <- ggplot(df) +
        geom_boxplot(aes(x = group, y = prop, color = group, fill = group), width = 0.9, position = position_dodge(width = 1), outlier.colour = NA) +
        geom_point(aes(x = group, y = prop, color = group, shape = data), size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 1)) +
        theme_bw() +
        ylab("Frequency") +
        xlab("") +
        theme(axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 14, face = "bold", color = "black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
          legend.title = element_blank(), legend.position = "right", 
          legend.key = element_blank(), legend.text = element_text(size = 14),
          strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0, face = "bold")) +
        scale_color_manual(values = color_groups) +
        scale_fill_manual(values = color_groupsb) +
        facet_wrap(~ cluster, scales = "free", ncol = nr_in_one_row)
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_", days[i] ,"2.pdf")), w = ww, h = hh)
      print(ggp)
      dev.off()
      
    }
    
    
    # ------------------------------------
    # plot for each sample as a bar with cluster composition
    
    days <- levels(ggdf$day)
    hh <- 5
    ww <- 10
    
    for(i in 1:nlevels(ggdf$day)){
      # i = 1
      
      df <- ggdf[ggdf$day == days[i], , drop = FALSE]
      
      ggp <- ggplot(df, aes(x = samp, y = prop, fill = cluster)) +
        geom_bar(stat = "identity") +
        theme_bw() +
        ylab("Frequency") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12),
          axis.title.y = element_text(size = 12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
          legend.title = element_blank(), legend.position = "right", legend.key = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0)) +
        scale_fill_manual(values = colors_clusters) +
        facet_wrap(~ group, scales = "free_x") 
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_", days[i] ,"3.pdf")), w = ww, h = hh)
      print(ggp)
      dev.off()
      
    }
    
    
  }
  
  
  
  
  return(NULL)
  
}



