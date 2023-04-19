# shru
SHared R Utilities for behavioural genomics and beyond

Includes the previous semPlate package utilities for generating path diagrams from fitted SEM model results.

## Included utilities
- Supermunge: A function to munge (QC, harmonise to reference) GWAS summary statistics or generic text based table-like genetic variant file data, like .bim files for example. Parses tab-separated files or inline data frames. Can impute missing variant information using the novel LD-IMP method. Can lift-over genetic coordinates to new coordinates by using a UCSC chain file.

- Semplate: Various functions for parsing and visualising SEM data, Genomic SEM in particular.

## Dependencies
You may have to load either tidyverse and/or DiagrammeR for the functionality that requires it.

## Code usage examples

Supermunge with basic options; perform munge using the provided reference variant list and apply basic QC
```[r]
munge$results <- supermunge(
                filePaths = munge$filesToUse, #list of file paths
                refFilePath = filepath.SNPReference.1kg, #reference variant list (SNP,CHR,BP, A1, A2, MAF, optional CM and L2)
                traitNames = munge$traitNamesToUse, #names/codes of the traits in the filePaths
                N = munge$NToUse, #N values per train in filePaths to use in case the sumstats do not contain N
                pathDirOutput = folderpath.data.sumstats.munged #path to the folder where to put the output files, named by the (trait name).gz
    )
```

Supermunge with no munge-step but file interpretation only and specific cleaning filters
```[r]
munge$results <- supermunge(
                filePaths = munge$filesToUse, #list of file paths
                traitNames = munge$traitNamesToUse, #names/codes of the traits in the filePaths
                N = munge$NToUse, #N values per train in filePaths to use in case the sumstats do not contain N
                pathDirOutput = folderpath.data.sumstats.munged, #path to the folder where to put the output files, named by the (trait name).gz
                process=F,
                filter.maf = 0.01,
                filter.info = 0.6,
                filter.chr = list(23,24,25,26) #remove chromosomes 23, 24, 25 and 26 using PLINK numeric format
    )
```

Supermunge with LD-IMP GWAS summary statistics imputation and other advanced features activated
```[r]
    munge$results <- supermunge(
                filePaths = munge$filesToUse, #list of file paths
                rsSynonymsFilePath = filepath.rsSynonyms.dbSNP151, #path to rs-id synonym file (each row has a list of id’s, separated by space, with the first id to be the preferred id)
                refFilePath = filepath.SNPReference.1kg, #reference variant list (SNP,CHR,BP, A1, A2, MAF, optional CM and L2)
                ldDirPath = folderpath.data.mvLDSC.ld.1kg, #if the variant list does not contain ld-scores (L2) you need to specify a folder to read LD-scores from – used for LD-IMP
                traitNames = munge$traitNamesToUse, #names/codes of the traits in the filePaths
                chainFilePath = file.path(p$folderpath.data,"alignment_chains","hg19ToHg38.over.chain.gz"), #optional: chain file to perform liftover if necessary
                N = munge$NToUse, #N values per train in filePaths to use in case the sumstats do not contain N
                imputeFromLD=T, #perform LD-IMP
                imputeFrameLenBp = 500000,
                imputeFrameLenCM=NULL, #the default setting is to use the 0.5 cM setting here. You will however need to have cM values in your variant list. Set this to NULL to use the setting in imputeFrameLenBp instead.
                pathDirOutput = folderpath.data.sumstats.munged #path to the folder where to put the output files, named by the (trait name).gz
    )

