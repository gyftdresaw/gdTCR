

# counting/identifying duplicated sequences -- without repeats or double counting
setwd('~/Documents/Dinner/gd TCR/TCR_analysis 062116/gdTCR/')

load('data/all_patients_split.RData')

## getting all duplicated CDR3 sequences ##
DT = all_patients_split.df[,c('Chain','CDR3')]
dup_cdr3_set = unique(DT[duplicated(DT),])
dup_cdr3 = list()
for (di in 1:nrow(dup_cdr3_set)) {
  dup_cdr3 = c(dup_cdr3,list(all_patients_split.df[with(all_patients_split.df,
                                                        Chain==dup_cdr3_set$Chain[di] & CDR3==dup_cdr3_set$CDR3[di]),]))
}


