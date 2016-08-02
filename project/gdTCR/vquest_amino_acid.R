
# parsing vquest by amino acid sequence matching
setwd('~/Documents/Dinner/gd TCR/TCR_analysis-master 062116/gdTCR/')

require(Biostrings) ## longest common prefix/suffix functions

# preloaded dataset
all_patients_pre.df = read.csv("data/data_summary.csv")

## 

vstring = ''
