# UntargetedMetabolomics_FoodGradeMicroorganisms_Fermentation

This repository contains the R-codes to perform the elaborations related to the work:

**Untargeted Metabolomic Profiling Reveals Multifunctional Bioactive Postbiotics in Cell-Free Supernatants from Food-Grade Microbial Fermentation Waste.** 
Denise Drago, Luca Marini, Eleonora Bafile, Federica Federici, Radmila Pavlovic, Marynka Ulaszewska, Gianfranco Frigerio, Laura Manna, Annapaola Andolfo.



The input tables can be obtained from the Zenodo repository:

https://doi.org/10.5281/zenodo.19251430



In order to run the elaborations with the same R-package versions, run the following code:

```r
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::restore()
```