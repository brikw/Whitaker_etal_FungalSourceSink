# "Local plants, not soils, are the primary source of foliar fungal community assembly in a C4 grass"
### Whitaker, B.K., Giauque, H., Timmerman, C., Birk, N., Hawkes, C.V.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5076546.svg)](https://doi.org/10.5281/zenodo.5076546)

```
Briana K. Whitaker. (2021, July 6). brikw/Whitaker_etal_FungalSourceSink: Local plants, 
not soils, are the primary source of foliar fungal community assembly in a C4 grass 
(Version 1.3). Zenodo. https://doi.org/10.5281/zenodo.5076546
```

This repository includes the R code, data files, small scripts, and a metadata file to supplement the manuscript by Whitaker et al. "Local plants, not soils, are the primary source of foliar fungal community assembly in a C4 grass". 

The sample matrix ("TXEG_Metadata_Combined.csv") includes information about each sample collected from
across a steep precipitaton gradient in central Texas. Sample types included a focal host at every site, Panicum hallii, each of a locally abundant C4 grass and dicot, as well as soil. Associated environmental measurements are included for the final sample set included in the manuscript.

The ASV matrix ("TXEG_SbyS_filtered.csv") includes the number of amplicon sequence variant (ASV) reads for each sample (rows) and unique ASV (columns). The raw sequence data used to generate this matrix is available through the NCBI SRA database under accession number PRJNA672680 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA672680/).

The taxonomy matrix ("TXEG_filteredSeqs_Classified.csv") includes the taxonomic classification of each unique ASV as was used in the manuscript.

Detailed information about column headers in the sample matrix and taxonomy matrix can be found in the metadata file "TXEG_Meta.xlsx".

THe three Rmd documents provide the pipeline for bioinformatic analysis using DADA2, taxonomic classification of fungal ASVs and removal of plant ASVs, and lastly the analysis of the fungal community as it appears in the paper. They should be run in this order: DADA2_pipeline.Rmd > TaxonomicClassification.Rmd > CommunityAnalyses.Rmd 

The /code folder contains the batch and R scripts necessary to run the full sourctracker analysis on a supercomputer. Data necessary to run the analyses can be found in /data folder.


Please see the manuscript for details and full reference information.

