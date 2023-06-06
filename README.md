# shru
SHared R Utilities for behavioural genomics and beyond

Includes the previous semPlate package utilities for generating path diagrams from fitted SEM model results.

## Included utilities
- Supermunge: A function to munge (QC, harmonise to reference) GWAS summary statistics or generic text based table-like genetic variant file data, like .bim files for example. Parses tab-separated files or inline data frames. Can impute missing variant information using the novel LD-IMP method. Can lift-over genetic coordinates to new coordinates by using a UCSC chain file.

- Semplate: Various functions for parsing and visualising SEM data, Genomic SEM in particular.

## Dependencies
You may have to additionally load either tidyverse and/or DiagrammeR for the functionality that requires it. These were left out of the package dependencies to increase the compatibility of SHRU with computational HPC environments that may not have the installed software to support some graphical routines.
