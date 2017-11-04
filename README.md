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

Analysis is done using R version 3.3.0 and Bioconductor version 3.3

To generate t-SNE plots with cells detected by CellCnn run:

```
make -f 000_cellcnn_plot_tsne_pipeline.mk
```



