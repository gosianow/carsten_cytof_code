

Sys.time()

# Load packages
library(gdata)
library(ggplot2)
library(reshape2)
library(limma) # for strsplit2
library(gtools) # for logit
library(plyr) # for rbind.fill
library(tools)
library(GGally)
library(ComplexHeatmap)

##############################################################################
# Test arguments
##############################################################################


outdir='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/correlation'
path_metadata=c('../carsten_cytof/PD1_project/CK_metadata/metadata_23_03.xlsx','../carsten_cytof/PD1_project/CK_metadata/metadata_29_03.xlsx')

path_frequencies=c('../carsten_cytof/PD1_project/CK_2016-06-23_03/050_frequencies/23_03_pca1_merging5_frequencies.xls','../carsten_cytof/PD1_project/CK_2016-06-29_03/050_frequencies/29_03_pca1_merging2_frequencies.xls')

## CD4
prefix='02_IFN_myeloid_CD4_'
path_bimatrix_frequencies='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD4/090_cytokine_bimatrix/23CD4allall_29CD4allall_02CD4v2_frequencies.xls'

## CD8
# prefix='02_IFN_myeloid_CD8_'
# path_bimatrix_frequencies='../carsten_cytof/PD1_project/CK_2016-06-merged_23_29/02v2_CD8/090_cytokine_bimatrix/23CD8allall_29CD8allall_02CD8v2_frequencies.xls'


##############################################################################
# Read in the arguments
##############################################################################

# rm(list = ls())
# 
# args <- (commandArgs(trailingOnly = TRUE))
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# 
# cat(paste0(args, collapse = "\n"), fill = TRUE)


##############################################################################

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- lapply(1:length(path_metadata), function(i){
  md <- read.xls(path_metadata[i], stringsAsFactors = FALSE)
  md
})

md <- plyr::rbind.fill(md)
rownames(md) <- md$shortname

### Factor arrangment
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)
md$data <- factor(md$data)
md$data_day <- interaction(md$data, md$day, lex.order = TRUE, drop = TRUE)

### Colors 
colors <- unique(md[, c("condition", "color")])
colors$condition <- factor(colors$condition)
## replace _ with \n
levels(colors$condition) <- gsub("_", "\n", levels(colors$condition ))

color_groups <- colors$color
names(color_groups) <- colors$condition

color_groupsb <- adjustcolor(color_groups, alpha = 0.3)
names(color_groupsb) <- colors$condition

color_samples <- md$color
names(color_samples) <- md$shortname

colors <- unique(md[, c("response", "color")])
color_response <- colors$color
names(color_response) <- colors$response


# ------------------------------------------------------------
# Load cluster frequencies
# ------------------------------------------------------------

prop <- lapply(1:length(path_frequencies), function(i){
  # i = 1
  prop <- read.table(path_frequencies[i], header = TRUE, sep = "\t", as.is = TRUE)
  print(prop[, c("cluster", "label")])
  return(prop)
})


prop_out <- Reduce(function(...) merge(..., by = c("cluster", "label"), all=TRUE, sort = FALSE), prop)
prop_out[, c("cluster", "label")]

## drop the 'drop' cluster
prop_out <- prop_out[prop_out$label != "drop", , drop = FALSE]

if(!all(complete.cases(prop_out))){
  stop("There are some clusters that are not common in the merged data sets or have different cluster number!")
}

## keep only those samples that are also in the metadata file
prop_out <- prop_out[, c("cluster", "label", md$shortname)]

prop_out <- prop_out[order(prop_out$cluster), , drop = FALSE]


### Create labels 
labels <- data.frame(cluster = prop_out$cluster, label = prop_out$label)

labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
labels


prop_tmp <- t(prop_out[, md$shortname])
colnames(prop_tmp) <- prop_out$label

prop <- data.frame(sample = md$shortname, prop_tmp, check.names = FALSE)

colnames(prop) <- gsub("+", "pos", colnames(prop), fixed = TRUE)
colnames(prop) <- gsub("-", "_", colnames(prop), fixed = TRUE)


# ------------------------------------------------------------
# Load the bimatrix frequencies
# ------------------------------------------------------------


biprop <- lapply(1:length(path_bimatrix_frequencies), function(i){
  # i = 1
  biprop <- read.table(path_bimatrix_frequencies[i], header = TRUE, sep = "\t", as.is = TRUE)
})

biprop_out <- Reduce(function(...) merge(..., by = c("cluster", "label"), all=TRUE, sort = FALSE), biprop)

biprop_tmp <- t(biprop_out[, md$shortname])
colnames(biprop_tmp) <- biprop_out$label

biprop <- data.frame(sample = md$shortname, biprop_tmp, check.names = FALSE)

colnames(biprop) <- gsub("+", "pos", colnames(biprop), fixed = TRUE)
colnames(biprop) <- gsub("-", "_", colnames(biprop), fixed = TRUE)

# -----------------------------------------------------------------------------
# Prepare the ggadf and gglabels objects
# -----------------------------------------------------------------------------

gglabels <- c("CD14pos_monos", "IFN_gpos")

### Prepare the data for plotting with ggplot
ggadf <- merge(biprop, prop, by = "sample", sort = FALSE)
rownames(ggadf) <- ggadf$sample


ggadf$group <- factor(md[rownames(ggadf), "condition"])
## replace _ with \n
levels(ggadf$group) <- gsub("_", "\n", levels(ggadf$group))
ggadf$data_day <- md[rownames(ggadf), "data_day"]
ggadf$data <- md[rownames(ggadf), "data"]

head(ggadf)


# -----------------------------------------------------------------------------
# Correlation analysis
# -----------------------------------------------------------------------------



### Pairs plot with GGally

# ggp <- ggpairs(ggadf[, gglabels]) +
#   theme_bw()
# 
# pdf(file.path(outdir, paste0(prefix, "frequencies_plot_corr_ggally.pdf")), 10, 10)
# print(ggp)
# dev.off()


shape_data_day <- c(19, 1, 17, 2)
names(shape_data_day) <- levels(ggadf$data_day)
shape_data_day


### Individual paired plots 

for(i in 1:(length(gglabels) - 1)){
  
  for(j in (i+1):length(gglabels)){
    # i = 1; j = 2
    
    ggp <- ggplot(ggadf) +
      geom_point(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]), shape = "data_day", color = "group"), size = 3, alpha = 0.8) +
      geom_smooth(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]))) +
      theme_bw() +
      theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.title = element_blank()) +
      scale_shape_manual(values = shape_data_day) +
      scale_color_manual(values = color_groups)
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot_corr_pairs_", gglabels[i], "_", gglabels[j] ,"_dataALL.pdf")), width = 7, height = 5)
    print(ggp)
    dev.off()
    
    ### Plots per data day
    for(dd in levels(ggadf$data_day)){
      
      ggp <- ggplot(ggadf[ggadf$data_day == dd, , drop = FALSE]) +
        geom_point(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]), shape = "data_day", color = "group"), size = 3, alpha = 0.8) +
        geom_smooth(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]))) +
        theme_bw() +
        theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 16), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          legend.title = element_blank()) +
        scale_shape_manual(values = shape_data_day) +
        scale_color_manual(values = color_groups)
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_corr_pairs_", gglabels[i], "_", gglabels[j], "_", gsub(".", "", dd, fixed = TRUE) ,".pdf")), width = 7, height = 5)
      print(ggp)
      dev.off()
      
    }
    
    
  }
  
}




### Heatmap with correlation 

corr_methods <- c("spearman")

for(m in 1:length(corr_methods)){
  # m = 1
  
  mat <- cor(ggadf[, gglabels], method = corr_methods[m], use = "complete.obs")
  
  out <- data.frame(cluster = rownames(mat), mat, stringsAsFactors = FALSE)
  
  write.table(out, file.path(outdir, paste0(prefix, "frequencies_plot_corr_corr_", corr_methods[m], "_dataALL.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  mat[upper.tri(mat)] <- NA
  diag(mat) <- NA
  
  if(length(gglabels) > 2){
    ### Using ComplexHeatmap
    legend_breaks = seq(from = -round(1), to = round(1), by = 0.5)
    
    ht1 <- Heatmap(mat, name = "Correlation", col = colorRampPalette(c("#dc143c", "#f5f5f5", "#4682b4"), space = "Lab")(15), na_col = "white", cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left", heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous", legend_height = unit(40, "mm")), cell_fun = function(j, i, x, y, w, h, col){
      if(j < i)
        grid.text(round(mat[i, j], 2), x, y)
    }) 
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot_corr_heat_", corr_methods[m] ,"_asis", "_dataALL.pdf")))
    draw(ht1)
    dev.off()
  }
  
  
  ### Analysis per data day
  for(dd in levels(ggadf$data_day)){
    
    mat <- cor(ggadf[ggadf$data_day == dd, gglabels, drop = FALSE], method = corr_methods[m], use = "complete.obs")
    
    out <- data.frame(cluster = rownames(mat), mat, stringsAsFactors = FALSE)
    
    write.table(out, file.path(outdir, paste0(prefix, "frequencies_plot_corr_corr_", corr_methods[m], "_", gsub(".", "", dd, fixed = TRUE) ,".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    
  }
  
  
  
}




### Test the correlation coefficients


corr_methods <- c("spearman")

for(m in 1:length(corr_methods)){
  # m = 1
  
  corr_pvs <- matrix(NA, length(gglabels), length(gglabels))
  colnames(corr_pvs) <- gglabels
  rownames(corr_pvs) <- gglabels
  
  
  for(i in 1:(length(gglabels) - 1)){
    
    for(j in (i+1):length(gglabels)){
      # i = 1; j = 2
      print(paste0(gglabels[i], " vs ", gglabels[j]))
      
      out <- cor.test(x = ggadf[, gglabels[i]], y = ggadf[, gglabels[j]], alternative = "two.sided", method = corr_methods[m])
      corr_pvs[j, i] <- out$p.value
      
      ggp_title <- paste0(round(out$estimate, 2))
      
      ### Plot
      ggp <- ggplot(ggadf) +
        geom_point(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]), shape = "data_day", color = "group"), size = 3, alpha = 0.8) +
        geom_smooth(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]))) +
        ggtitle(ggp_title) +
        theme_bw() +
        theme(title = element_text(size = 16), 
          axis.text = element_text(size = 14), 
          axis.title = element_text(size = 16), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          legend.title = element_blank()) +
        scale_shape_manual(values = shape_data_day) +
        scale_color_manual(values = color_groups)
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_plot_corr_pairs_", corr_methods[m], "_", gglabels[i], "_", gglabels[j] ,"_dataALL.pdf")), width = 7, height = 5)
      print(ggp)
      dev.off()
      
      
    }
    
  }
  
  corr_apvs <- matrix(p.adjust(corr_pvs, method = "BH"), nrow = length(gglabels), byrow = FALSE)
  colnames(corr_apvs) <- gglabels
  rownames(corr_apvs) <- gglabels
  
  out <- data.frame(cluster = rownames(corr_apvs), corr_apvs, stringsAsFactors = FALSE)
  
  write.table(out, file.path(outdir, paste0(prefix, "frequencies_plot_corr_apvs_", corr_methods[m], "_dataALL.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  out <- data.frame(cluster = rownames(corr_pvs), corr_pvs, stringsAsFactors = FALSE)
  
  write.table(out, file.path(outdir, paste0(prefix, "frequencies_plot_corr_pvs_", corr_methods[m], "_dataALL.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  
  if(length(gglabels) > 2){
    ### Using ComplexHeatmap
    legend_breaks = c(0, 0.05, 0.1, 1)
    
    ht1 <- Heatmap(corr_apvs, name = "Correlation", col = colorRampPalette(c("grey50", "grey70", "grey90", "grey90"))(4), na_col = "white", cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left", heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "discrete", legend_height = unit(40, "mm")), cell_fun = function(j, i, x, y, w, h, col){
      if(j < i)
        grid.text(sprintf("%.02e", corr_apvs[i, j]), x, y)
    }) 
    
    pdf(file.path(outdir, paste0(prefix, "frequencies_plot_corr_apvs_", corr_methods[m] ,"_asis", "_dataALL.pdf")))
    draw(ht1)
    dev.off()
  }
  
  
  ### Analysis per data day
  for(dd in levels(ggadf$data_day)){
    
    corr_pvs <- matrix(NA, length(gglabels), length(gglabels))
    colnames(corr_pvs) <- gglabels
    rownames(corr_pvs) <- gglabels
    
    for(i in 1:(length(gglabels) - 1)){
      
      for(j in (i+1):length(gglabels)){
        # i = 1; j = 2
        print(paste0(gglabels[i], " vs ", gglabels[j]))
        
        out <- cor.test(x = ggadf[ggadf$data_day == dd, gglabels[i]], y = ggadf[ggadf$data_day == dd, gglabels[j]], alternative = "two.sided", method = corr_methods[m])
        corr_pvs[j, i] <- out$p.value
        
        
        ggp_title <- paste0(round(out$estimate, 2))
        
        ### Plot 
        ggp <- ggplot(ggadf[ggadf$data_day == dd, , drop = FALSE]) +
          geom_point(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]), shape = "data_day", color = "group"), size = 3, alpha = 0.8) +
          geom_smooth(aes_string(x = as.character(gglabels[i]), y = as.character(gglabels[j]))) +
          ggtitle(ggp_title) +
          theme_bw() +
          theme(title = element_text(size = 16),
            axis.text = element_text(size = 14), 
            axis.title = element_text(size = 16), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(), 
            axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), 
            axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
            legend.title = element_blank()) +
          scale_shape_manual(values = shape_data_day) +
          scale_color_manual(values = color_groups)
        
        pdf(file.path(outdir, paste0(prefix, "frequencies_plot_corr_pairs_", corr_methods[m], "_", gglabels[i], "_", gglabels[j], "_", gsub(".", "", dd, fixed = TRUE) ,".pdf")), width = 7, height = 5)
        print(ggp)
        dev.off()
        
        
      }
      
    }
    
    corr_apvs <- matrix(p.adjust(corr_pvs, method = "BH"), nrow = length(gglabels), byrow = FALSE)
    colnames(corr_apvs) <- gglabels
    rownames(corr_apvs) <- gglabels
    
    out <- data.frame(cluster = rownames(corr_apvs), corr_apvs, stringsAsFactors = FALSE)
    
    write.table(out, file.path(outdir, paste0(prefix, "frequencies_plot_corr_apvs_", corr_methods[m], "_", gsub(".", "", dd, fixed = TRUE), ".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    out <- data.frame(cluster = rownames(corr_pvs), corr_pvs, stringsAsFactors = FALSE)
    
    write.table(out, file.path(outdir, paste0(prefix, "frequencies_plot_corr_pvs_", corr_methods[m], "_", gsub(".", "", dd, fixed = TRUE), ".xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    
  }
  
  
  
  
  
}




sessionInfo()





