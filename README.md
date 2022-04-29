# esMethods
This repository contains the data files and codes used in the paper:

Sasmita, K., & Swallow, K. M.(2022). Measuring event segmentation: An investigation into the stability of event boundary agreement across groups. *Behavior Research Methods* doi: https://doi.org/10.3758/s13428-022-01832-5

# Directory 

Folder: data/bootstrapped
>Contains all the bootstrapped agreement measures organized by each measure 
  
Folder: data/raw 
>Contains formatted participant segmentation data
>>'*_SegmentData_Filtered.txt' - segmentation data organized as each participants' button press times  
>>'*_TimeSeries_Filtered.txt' -  segmentation data organized as time series of presence or absence of button press every 1s bin

Folder: code
>Contains all the codes used to generate bootstraps and analyses
>>'get.esMethods.metrices.R' - contain all the functions to calculate the various agreement measures, called by 'esMethods_GenerateIterations' and 'esMethods_getFullSampleAgreement.rmd'
>>
>>'esMethods_GenerateIterations.rmd' - generate bootstraps for all agreement measures for both datasets, generated output goes to ../data/bootstrapped
>>
>>'esMethods_ManuscriptAnalysis.rmd' - perform all analyses and plotting included in the manuscript, plot output goes to ../plots
>>'esMethods_getFullSampleAgreement.rmd' - calculate the agreement values for the full samples (commercial-lab and everyday-online) used for analyses. This script can act as an example on how to use the functions in get.esMethods.metrics.R
>>
>>'esMethods_getBootstrapStatistics.rmd' - calculate the descriptive statistics (mean of sample means, sd of sample means, sd of population mean, 2.5 & 97.5 percentile of sample means) of bootstrapped agreement values for all sample sizes (n = 2 to 32). Output: .csv files for each agreement measure in each dataset, generated output goes to ../data/summary_statistics

# Run code 

To run all data processing and analyses reported in the manuscript: 
1. Generate bootstrapped data: run 'esMethods_GenerateIterations.rmd'
2. Run analyses and generate plots: run 'esMethods_ManuscriptAnalysis.rmd'
3. Calculate agreement values for the full samples: run 'esMethods_getFullSampleAgreement.rmd'
4. Calculate descriptive statistics for the bootstrapped agreement values: run 'esMethods_getBootstrapStatistics.rmd'
