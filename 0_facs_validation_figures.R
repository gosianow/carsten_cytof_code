library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(limma) # for strsplit2
library(RColorBrewer)
library(pheatmap)
library(gtools) # for logit

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/FACS_validation'
path_data <- "FACSvalidation_01_GN.xlsx"

setwd(rwd)

prefix <- "facs_validation_"
suffix <- ""
outdir <- "./"

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


### Colors
color_groups <- c("#000000", "#CC79A7", "#009E73")
names(color_groups) <- c("HD", "NR", "R")

color_groupsb <- adjustcolor(color_groups, alpha = 0.3)
names(color_groupsb) <- c("HD", "NR", "R")

### Read in data for plotting
data <- read.xls(path_data, stringsAsFactors=FALSE)

ggdf <- data

ggdf$group <- factor(ggdf$response, levels = c("HD", "NR", "R"))
ggdf$cluster <- factor(ggdf$cell_type)
ggdf$prop <- ggdf$frequency



# ------------------------------------
# plot all clusters in one pdf; colors per group; boxplots + points


ggp <- ggplot(ggdf, aes(x = cluster, y = prop, color = group, fill = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(size = 2, shape = 16, alpha = 0.8, position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.7)) +
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
  scale_fill_manual(values = color_groupsb)

pdf(file.path(outdir, paste0(prefix, "frequencies_plot.pdf")), w = nlevels(ggdf$cluster) + 1, h=4)
print(ggp)
dev.off()



