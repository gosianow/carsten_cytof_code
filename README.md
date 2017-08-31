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
