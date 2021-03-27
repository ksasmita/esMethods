# esMethods
This repository contains the data files and codes used in the paper:

Measuring Events: An Investigation into the Stability of Event Segmentation Pattens Across Groups
by K.Sasmita and K.M.Swallow 2021

# Directory 

Folder: data/bootstrapped
>Contains all the bootstrapped agreement measures organized by each measure 
  
Folder: data/raw 
>Contains formatted participant segmentation data
>>'*_SegmentData.txt' - segmentation data organized as each participants' button press times  
>>'*_TimeSeries.txt' -  segmentation data organized as time series of presence or absence of button press every 1s bin

Folder: code
>Contains all the codes used to generate bootstraps and analyses
>>'get.esMethods.metrices.R' - contain all the functions to calculate the various agreement measures, called by 'esMethods_GenerateIterations'
>>
>>'esMethods_GenerateIterations.rmd' - generated bootstraps for all agreement measures for both datasets, generated output goes to ../data/bootstrapped
>>
>>'esMethods_ManuscriptAnalysis.rmd' - perform all analyses and plotting included in the manuscript, plot output goes to ../plots

# Run code 

To run all data processing and analyses reported in the manuscript: 
1. Generate bootstrapped data: run 'esMethods_GenerateIterations.rmd'
2. Run analyses: run 'esMethods_ManuscriptAnalysis.rmd'
