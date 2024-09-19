RVTN Within Host Evolution of SARS-CoV-2
--
This repository contains code and small intermediate data files associated with the manuscript "In depth sequencing of a serially sampled household cohort reveals the within-host dynamics of Omicron SARS-CoV-2 and rare selection of novel spike variants". 

The directory is organized as follows:

**Pipeline** -  This folder contains shell scripts for processing raw fastq files and calling variants. See manuscript for details of sample processing.  Raw sequencing data is available in the SRA as BioProject PRJNA1159790.  

**References** - This folder contains the reference sequence, gff file, and information needed for divergence rate estimates. 

**Metadata** - This folder contains metadata for all specimens. Meta_1000.csv contains the meta data for the samples with high quality sequencing used for analysis. 

**Scripts** - This folder contains custom R scripts to perform analyses in the manuscript. The files are numbered in order of usage with figure scripts being run last. Most intermediate files are included in this repo and each script should be able to run independently of each other. 

**Results** - This file contains the variant file generated from ivar with associated metadata. The variants have been fully filtered and is the input file for all downstram analyses
This folder also contains intermediate files and results of analyses. 
