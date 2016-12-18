


plot_frequencies <- function(ggdf, color_groups, color_groupsb, outdir, prefix, pdf_hight){
  
  
  if(!all(is.na(ggdf$prop))){
    
    # ------------------------------------
    # plot all clusters in one pdf; colors per response; boxplots + points
    
    
    ggp <- ggplot(ggdf, aes(x = cluster, y = prop, color = group, shape = data, fill = group)) +
      geom_boxplot(width = 0.95, position = position_dodge(width = 0.95), outlier.colour = NA) +
      geom_point(size=1.5, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
      theme_bw() +
      ylab("Frequency") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12, face="bold"), 
        axis.title.y = element_text(size=12, face="bold"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
        axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
        legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = color_groups) +
      scale_fill_manual(values = color_groupsb) +
      facet_wrap(~ day)
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")), w = nlevels(ggdf$cluster) + 3, h = pdf_hight)
    print(ggp)
    dev.off()
    
    
    # ------------------------------------
    ### plot each cluster as a separate page in the pdf file
    ggp <- list()
    clusters <- levels(ggdf$cluster)
    
    for(i in 1:nlevels(ggdf$cluster)){
      # i = 1
      
      df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
      
      ggp[[i]] <- ggplot(df, aes(x = group, y = prop, color = group, shape = data, fill = group)) +
        geom_boxplot(width = 0.95, position = position_dodge(width = 0.95), outlier.colour = NA) +
        geom_point(size=3, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
        ggtitle(clusters[i]) +
        theme_bw() +
        ylab("Frequency") +
        xlab("") +
        theme(axis.text.x = element_text(size=12, face="bold"), 
          axis.title.y = element_text(size=12, face="bold"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          legend.position = "none") +
        guides(color = guide_legend(ncol = 1)) +
        scale_color_manual(values = color_groups) +
        scale_fill_manual(values = color_groupsb)
      
    }
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot2.pdf")), w=5, h=4, onefile=TRUE)
    for(i in seq(length(ggp)))
      print(ggp[[i]])
    dev.off()
    
    
    # ------------------------------------
    # plot all clusters in one pdf; colors per response; separate pdf for base and tx; one boxplot + points
    
    
    days <- levels(ggdf$day)
    
    for(i in 1:nlevels(ggdf$day)){
      # i = 1
      
      df <- ggdf[ggdf$day == days[i], , drop = FALSE]
      
      ggp <- ggplot(df) +
        geom_boxplot(aes(x = cluster, y = prop, color = group, fill = group), width = 0.95, position = position_dodge(width = 0.95), outlier.colour = NA) +
        geom_point(aes(x = cluster, y = prop, color = group, shape = data), size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
        theme_bw() +
        ylab("Frequency") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12, face="bold"), 
          axis.title.y = element_text(size=12, face="bold"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
          legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
        scale_color_manual(values = color_groups) +
        scale_fill_manual(values = color_groupsb)
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_boxplotpoints_", days[i] ,".pdf")), w = nlevels(ggdf$cluster)/2 + 3, h = pdf_hight)
      print(ggp)
      dev.off()
      
    }
    
    
    # ------------------------------------
    # plot all clusters in one pdf; colors per response; separate pdf for base and tx; one boxplot + points + facet per cluster
    
    days <- levels(ggdf$day)
    
    nr_cluster <- nlevels(ggdf$cluster)
    nrow <- ifelse(nr_cluster < 6, 1, 2)
    h <- ifelse(nr_cluster < 6, 2, 4)
    w <- ifelse(nr_cluster < 6, nr_cluster * 1.25 + 2, ceiling(nr_cluster/2) * 1.5 + 2)
    
    for(i in 1:nlevels(ggdf$day)){
      # i = 1
      
      df <- ggdf[ggdf$day == days[i], , drop = FALSE]
      
      ggp <- ggplot(df) +
        geom_boxplot(aes(x = cluster, y = prop, color = group, fill = group), width = 0.95, position = position_dodge(width = 0.95), outlier.colour = NA) +
        geom_point(aes(x = cluster, y = prop, color = group, shape = data), size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
        theme_bw() +
        ylab("Frequency") +
        xlab("") +
        theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, face="bold"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
          legend.title = element_blank(), legend.position = "right", legend.key = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size=9, hjust = 0)) +
        scale_color_manual(values = color_groups) +
        scale_fill_manual(values = color_groupsb) +
        facet_wrap(~ cluster, scales = "free", nrow = nrow)
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_boxplotpoints_", days[i] ,"2.pdf")), w = w, h = h)
      print(ggp)
      dev.off()
      
    }
    
    # ------------------------------------
    # plot all clusters in one pdf; colors per response; separate pdf for base and tx; one boxplot + points + pdf page per cluster
    
    days <- levels(ggdf$day)
    clusters <- levels(ggdf$cluster)
    
    
    for(i in 1:nlevels(ggdf$day)){
      # i = 1
      
      ggp <- list()
      
      for(j in 1:nlevels(ggdf$cluster)){
        
        df <- ggdf[ggdf$day == days[i] & ggdf$cluster == clusters[j], , drop = FALSE]
        
        ggp[[j]] <- ggplot(df) +
          geom_boxplot(aes(x = cluster, y = prop, color = group, fill = group), width = 0.95, position = position_dodge(width = 0.95), outlier.colour = NA) +
          geom_point(aes(x = cluster, y = prop, color = group, shape = data), size=2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.95, jitter.height = 0, dodge.width = 0.95)) +
          theme_bw() +
          ggtitle(clusters[j]) +
          ylab("Frequency") +
          xlab("") +
          theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12, face="bold"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(), 
            axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
            axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
            legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
          scale_color_manual(values = color_groups) +
          scale_fill_manual(values = color_groupsb)
      }
      
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_boxplotpoints_", days[i] ,"3.pdf")), w = 2.5, h = 3)
      for(j in seq(length(ggp)))
        print(ggp[[j]])
      dev.off()
      
    }
    
    
  }else{
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")))
    plot(1, type="n", axes=F, xlab="", ylab="")
    dev.off()
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot2.pdf")))
    plot(1, type="n", axes=F, xlab="", ylab="")
    dev.off()
    
    for(i in 1:nlevels(ggdf$day)){
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_boxplotpoints_", days[i] ,".pdf")))
      plot(1, type="n", axes=F, xlab="", ylab="")
      dev.off()
    }
    
    for(i in 1:nlevels(ggdf$day)){
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_boxplotpoints_", days[i] ,"2.pdf")))
      plot(1, type="n", axes=F, xlab="", ylab="")
      dev.off()
    }
    
    for(i in 1:nlevels(ggdf$day)){
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_boxplotpoints_", days[i] ,"3.pdf")))
      plot(1, type="n", axes=F, xlab="", ylab="")
      dev.off()
    }
    
  }
  
  
  
  return(NULL)
  
}



