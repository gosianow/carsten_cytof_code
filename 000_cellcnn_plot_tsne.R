

Sys.time()

# Load packages
library(flowCore)
library(gdata)
library(Rtsne)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(coop) # cosine
library(limma) 
library(tools)

##############################################################################
# Test arguments
##############################################################################

args <- NULL

prefix='23_01_pca1_combined_'
outdir='../PD1_CellCnn_Lukas_June2017/cellcnn_tsne'
path_metadata='../carsten_cytof/PD1_project/CK_metadata/metadata_23_01.xlsx'
path_rtsne_out='../carsten_cytof/PD1_project/CK_2016-06-23_01/040_tsnemaps/23_01_pca1_rtsne_out.rds'
path_rtsne_data='../carsten_cytof/PD1_project/CK_2016-06-23_01/040_tsnemaps/23_01_pca1_rtsne_data.xls'
dir_fcs='../carsten_cytof/PD1_project/CK_2016-06-23_01/010_cleanfcs'
dir_cellcnn_files='../PD1_CellCnn_Lukas_June2017/panel1_base_combined/selected_cells'
day='base'



##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)

##############################################################################


if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


suffix <- ""
pdf_width=15
pdf_height=7


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- read.xls(path_metadata, stringsAsFactors=FALSE)

# add more info about samples
cond_split <- strsplit2(md$condition, "_")
colnames(cond_split) <- c("day", "response")

md[, c("day", "response")] <- cond_split
md$response <- factor(md$response, levels = c("NR", "R", "HD"))
md$response <- factor(md$response)
md$day <- factor(md$day, levels = c("base", "tx"))
md$day <- factor(md$day)
md$patient_id <- factor(md$patient_id)


# ------------------------------------------------------------
# Load fcs files
# ------------------------------------------------------------

### Need to load them to find out the IDs of cells used in t-SNE 


# define FCS file names
f <- file.path(dir_fcs, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)

## Create sample info
samp <- rep(names(fcs), sapply(fcs, nrow))
cell_ids = 1:length(samp)

cell_ids_split <- split(cell_ids, factor(samp, levels = md$shortname))

sapply(cell_ids_split, length)

# ------------------------------------------------------------
# Load cell cnn results
# ------------------------------------------------------------

samps2keep <- md$day == day & md$response %in% c("R", "NR")

fcnn <- file.path(dir_cellcnn_files, paste0(file_path_sans_ext(md$filename[samps2keep]), "_transf_selected_cells.csv"))
names(fcnn) <- md$shortname[samps2keep]
fcnn


cnn_cells <- lapply(fcnn, function(x){
  
  sel_cells <- read.csv(x, header = TRUE, sep = ",")
  
  out <- factor(sel_cells[, 2])
  
  return(out)
  
})

sapply(cnn_cells, length)




# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

rtsne_out <- readRDS(path_rtsne_out)

rtsne_data <- read.table(path_rtsne_data, header = TRUE, sep = "\t", as.is = TRUE)

head(rtsne_data)



# ------------------------------------------------------------
# Identify cells to plot
# ------------------------------------------------------------

tsne_cells <- lapply(cell_ids_split, function(x){
  
  x %in% rtsne_data$cell_index
  
})


rtsne <- data.frame(tSNE1 = rtsne_out$Y[, 1], tSNE2 = rtsne_out$Y[, 2], sample = rtsne_data$sample_name)

rtsne_split <- split(rtsne, factor(rtsne_data$sample_name, levels = md$shortname))




ggdf <- lapply(names(cnn_cells), function(i){
  
  cnn_filter <- cnn_cells[[i]][tsne_cells[[i]]]
  
  out <- data.frame(rtsne_split[[i]], cnn_filter = cnn_filter)
  
  return(out)
})


ggdf <- plyr::rbind.fill(ggdf)


# add group info
mm <- match(ggdf$sample, md$shortname)
ggdf$group <- md$response[mm]


colors_filter <- c("#b2beb5", "#ff2052", "#9932cc", "#006a4e", "#ff8c00")[1:nlevels(ggdf$cnn_filter)]
names(colors_filter) <- levels(ggdf$cnn_filter)



# ----------------------------------------------------------------------
# tSNE plots
# ----------------------------------------------------------------------






# ----------------------------------------------------------------------
### Plot of tsne - all cells, all clusters



## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cnn_filter)) +
  geom_point(size = 1) +
  labs(x = "t-SNE1", y="t-SNE2")+ 
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
    axis.title  = element_text(size = 15),
    strip.text = element_text(size = 15, hjust = 0),
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_color_manual(values = colors_filter) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()




## facet per group
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cnn_filter)) +
  geom_point(size=1) +
  facet_wrap(~ group) +
  labs(x = "t-SNE1", y="t-SNE2")+
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
    axis.title  = element_text(size = 15),
    strip.text = element_text(size = 15, hjust = 0), 
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(size = 0.3, linetype = "solid", color = "black")) +
  scale_color_manual(values = colors_filter) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEgroup", suffix, ".pdf")), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



















sessionInfo()











