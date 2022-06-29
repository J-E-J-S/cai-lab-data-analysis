# PhD Data Analysis Scripts  
**Collection of scripts to analyse sequencing / machine data from the lab.**  
### Tools  
`conserved-genome-variance-bgi`: This script will analyse whole genome sequence data from BGI for a collection of strains and return a dataframe with all genetic variance with >1 occurences filtered to parental strain variance and for quality reads. This script will also output a dataframe for each strain with collated variance and filtered for parental strain and quality reads. Create a new variable with a path to each strains /SNV_Indel dir  
`plate-reader-biotek`: Create ribboned growth curves (OD600) normalized to background with Biotek plate reader.  
`plate-reader-tecan`: Create ribboned growth curves (OD600) normalized to background with Tecan plate reader.  
`FP-assay-biotek`: Create end-point GFP barcharts from the Biotek plate reader.
