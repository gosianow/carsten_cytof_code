


plot_heatmaps_for_sign_freqs <- function(expr_all, md, FDR_cutoff, pval_name2 = NULL, adjpval_name2 = NULL, pval_name_list, adjpval_name_list, breaks, legend_breaks, outdir, prefix, prefix2, suffix, color_response){
  # pval_name2, adjpval_name2 - p-value for plotting pheatmap2 - adjp_NRvsR contrast
  # pval_name_list, adjpval_name_list - p-value for plotting the pheatmap3 with everything
  
  rownames(md) <- md$shortname
  
  # -----------------------------
  ### Plot one heatmap with R vs NR
  ### pheatmap2
  
  ## order by p-value
  if(is.null(adjpval_name2)){
      which_top_pvs <- FALSE
    }else{
      adjpval_name <- adjpval_name2
      pval_name <- pval_name2
      expr_all <- expr_all[order(expr_all[, pval_name]), , drop = FALSE]
      
      which_top_pvs <- expr_all[, adjpval_name] < FDR_cutoff & !is.na(expr_all[, adjpval_name])
    }

  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0) {
    print("Plot pheatmap2")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- NULL
    
    ## expression scaled by row
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")")
    labels_col <- colnames(expr)
    
    annotation_col <- data.frame(response = factor(md[colnames(expr), "response"]))
    rownames(annotation_col) <- colnames(expr)
    
    annotation_colors <- list(response = color_response[levels(annotation_col$response)])
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap2", suffix, ".pdf")))
    
    
  }else{
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap2", suffix, ".pdf")))
    plot(1, type="n", axes=F, xlab="", ylab="")
    dev.off()
    
  }
  
  
  # -----------------------------
  ### Plot two heatmaps with R vs NR for base and tx
  
  for(i in levels(md$day)){
    # i = "base"
    print(paste0("Plot pheatmap_", i))
    
    adjpval_name <- paste0("adjp_NRvsR_", i)
    pval_name <- paste0("pval_NRvsR_", i)
    
    ## group the expression by cluster
    expr_all <- expr_all[order(expr_all[, pval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < FDR_cutoff & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
    
    if(sum(which_top_pvs) > 0){
      
      expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
      
      # -----------------------------
      ## order the samples by NR and R
      
      samples2plot <- md[md$response %in% c("NR", "R"), ]
      samples2plot <- samples2plot[order(samples2plot$response, samples2plot$day), ]
      samples2plot <- samples2plot[grep(i, samples2plot$day), ]
      samples2plot <- samples2plot$shortname
      
      ## gap in the heatmap
      gaps_col <- sum(grepl("_NR", samples2plot))
      gaps_row <- NULL
      
      ## expression scaled by row
      expr <- expr_heat[ , samples2plot, drop = FALSE]
      
      labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")")
      labels_col <- colnames(expr)
      
      annotation_col <- data.frame(response = factor(md[colnames(expr), "response"]))
      rownames(annotation_col) <- colnames(expr)
      
      annotation_colors <- list(response = color_response[levels(annotation_col$response)])
      
      pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap_", i, suffix, ".pdf")))
      
      
    }else{
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap_", i, suffix, ".pdf")))
      plot(1, type="n", axes=F, xlab="", ylab="")
      dev.off()
      
    }
    
  }
  
  # -----------------------------
  ### Plot one heatmap with R vs NR + heatmap with p-values for NRvsR_base, NRvsR_tx and NRvsR_basevstx
  ### pheatmap3
  
  ## group the expression by cluster and order by adjpval
  for(i in length(pval_name_list):1){
    expr_all <- expr_all[order(expr_all[, pval_name_list[i]]), , drop = FALSE]
  }
  
  
  # which_top_pvs <- rowSums(expr_all[, adjpval_name_list] < FDR_cutoff, na.rm = TRUE) > 0 & rowSums(is.na(expr_all[, adjpval_name_list])) == 0
  # which_top_pvs <- rowSums(is.na(expr_all[, adjpval_name_list])) < length(adjpval_name_list)
  which_top_pvs <- rep(TRUE, nrow(expr_all))
  which(which_top_pvs)
  
  if(sum(which_top_pvs) > 0){
    print("Plot pheatmap3")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    
    # -----------------------------
    ## order the samples by base and tx
    samples2plot <- md[md$response %in% c("NR", "R"), ]
    samples2plot <- samples2plot$shortname[order(samples2plot$day, samples2plot$response)]
    
    ## gap in the heatmap
    gaps_col <- c(max(grep("base_NR", samples2plot)), rep(max(grep("base", samples2plot)), 2), max(grep("tx_NR", samples2plot)))
    gaps_row <- NULL
    
    ## expression
    expr <- expr_heat[ , samples2plot, drop = FALSE]
    
    labels_row <- paste0(expr_heat$label)
    labels_col <- colnames(expr)
    
    annotation_col <- data.frame(response = factor(md[colnames(expr), "response"]))
    rownames(annotation_col) <- colnames(expr)
    
    annotation_colors <- list(response = color_response[levels(annotation_col$response)])
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap3", suffix, ".pdf")))
    
    
    pvs_heat <- expr_heat[, adjpval_name_list, drop = FALSE]
    
    labels_col <- colnames(pvs_heat)
    gaps_col <- NULL
    
    
    pheatmap(pvs_heat, cellwidth = 60, cellheight = 24, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap3pvs", suffix, ".pdf")))
    
    
  }
  

  
  return(NULL)
  
}



