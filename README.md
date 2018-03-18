# Analysis of the CyTOF data for the PD-1 project

Update the `RWD_MAIN` variable, which defines a path to the project directory, in the `Makefile`.

To run the pipeline cd to the directory that contains the code from this repository and type:

```
make
```

To update the time stamps of the files without rerunning the analysis (to touch the files) type:

```
make MAKEARGS="-t"
```

To generate t-SNE plots with cells detected by CellCnn run:

```
make -f 000_cellcnn_plot_tsne_pipeline.mk
```

## Dependencies 

Analysis is done using R version 3.3.0 with Bioconductor version 3.3

List of needed packages:

```
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
library(ggraph)
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

```

