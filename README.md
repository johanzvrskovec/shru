# shru

SHared R Utilities for behavioural genomics and beyond

## Included utilities

-   Supermunge: A function to munge (QC, harmonise to reference) GWAS summary statistics or generic text based table-like genetic variant file data, like .bim files for example. Parses tab-separated files or inline data frames.

    -   Can perform fixed-effect meta-analysis of multiple datasets using either inverse variance weighting or weighting by sample size.

    -   Can impute missing variant information using the novel LD-IMP method.

    -   Can lift-over genetic coordinates to new coordinates by using a UCSC chain file (still needs to be validated that it produces correct mappings - PROBABLY NOT FUNCTIONAL).

-   LDSC++: A modded version of the ldsc-function from the [Genomic SEM](https://github.com/GenomicSEM/GenomicSEM) package.

    -   Improved handling of GWAS summary statistics datasets with varying number of genetic variants and local genetic architecture, yielding more accurate estimates of genome-wide genetic covariance and its standard error.

    -   Compatible with Genomic SEM.

-   Semplate: Various functions for parsing and visualising SEM data, Genomic SEM in particular.

Please consider citing the pre-print here: <https://doi.org/10.1101/2025.04.07.25324946> if you use any of the routines within the package. The pre-print mainly describes the LDSC++ method, but also the routines in Supermunge.

## Dependencies

You may have to additionally load either tidyverse and/or DiagrammeR for the functionality that requires it. These were left out of the package dependencies to increase the compatibility of SHRU with computational HPC environments that may not have the installed software to support some graphical routines.
