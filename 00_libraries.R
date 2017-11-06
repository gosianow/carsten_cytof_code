### R version 3.3.0 (2016-05-03)
### Bioconductor version 3.3

### Various
library(gdata) # to read in xlsx read.xls()
library(limma) # for strsplit2
library(Repitools)
library(plyr)
library(coop)
library(gtools) # for logit
library(tools)

### Plotting
library(reshape2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(igraph)
library(GGally) # ggpair()
library(ggdendro)
library(UpSetR)
library(pheatmap)
library(ComplexHeatmap)

### Cytof methods
library(flowCore)
library(cytofkit)
library(citrus)
library(FlowSOM)
library(ConsensusClusterPlus)
library(Rtsne)

### Mixed models
library(lme4) # for fitting mixed models
library(multcomp) # for contrasts glht()

library(glmmADMB) # for glmmadmb() installed from http://glmmadmb.r-forge.r-project.org/
library(robustbase) # glmrob
library(robust) # glmRob
library(MASS) # rlm
library(betareg) # betareg() does not work with glht()


sessionInfo()

si <- sessionInfo()





