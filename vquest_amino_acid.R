
# parsing vquest by amino acid sequence matching
setwd('~/Documents/Dinner/gd TCR/TCR_analysis 062116/gdTCR/')

require(Biostrings) ## longest common prefix/suffix functions

# preloaded dataset
all_patients_pre.df = read.csv("data/data_summary.csv")

## v and j gene ends explicitly written out
v_genes.df = rbind(c('TRG','1','CATWDR'),
                   c('TRG','2','CATWDG'),
                   c('TRG','3','CATWDR'),
                   c('TRG','4','CATWDG'),
                   c('TRG','5','CATWDR'),
                   c('TRG','8','CATWDR'),
                   c('TRG','9','CALWEV'),
                   c('TRG','10','CAAWWV'),
                   c('TRD','1','CALGE'),
                   c('TRD','2','CACDT'),
                   c('TRD','3','CAF'),
                   c('TRD','TRAV9-2','CALS'),
                   c('TRD','TRAV21','CAVR'),
                   c('TRD','TRAV22','CAVE'),
                   c('TRD','TRAV41','CAVR'),
                   c('TRD','TRAV14/DV4','CAMRE'),
                   c('TRD','TRAV29/DV5','CAAS'),
                   c('TRD','TRAV38-2/DV8','CAYRS'))
v_genes.df = data.frame(v_genes.df,stringsAsFactors=F)
colnames(v_genes.df) = c('Chain','TRV','SEQ')

j_genes.df = rbind(c('TRG','1','NYYKKLF'),
                   c('TRG','1 or 2','NYYKKLF'),
                   c('TRG','2','NYYKKLF'),
                   c('TRG','2 or 1','NYYKKLF'),
                   c('TRG','P','GQELGKKIKVF'),
                   c('TRG','P or 1 or 2','GQELGKKIKVF'),
                   c('TRG','P1','TTGWFKIF'),
                   c('TRG','P2','SSDWIKTF'),
                   c('TRG','5','NYYKKLF'),
                   c('TRD','1','TDKLIF'),
                   c('TRD','2','LTAQLFF'),
                   c('TRD','3','SWDTRQMFF'),
                   c('TRD','4','RPLIF'))
j_genes.df = data.frame(j_genes.df,stringsAsFactors=F)
colnames(j_genes.df) = c('Chain','TRJ','SEQ')

split_cdr3 = function(x){
  vstring = v_genes.df$SEQ[which(v_genes.df$Chain == x['Chain'] & v_genes.df$TRV == x['TRV'])]
  jstring = j_genes.df$SEQ[which(j_genes.df$Chain == x['Chain'] & j_genes.df$TRJ == x['TRJ'])]
  cdr3string = x['CDR3']
  
  pre_index = lcprefix(vstring,cdr3string)
  suf_index = lcsuffix(jstring,cdr3string)

  split_items = c(substring(cdr3string,1,pre_index),
                  substring(cdr3string,pre_index+1,nchar(cdr3string)-suf_index),
                  substring(cdr3string,nchar(cdr3string)-suf_index+1))
  names(split_items) = c('CDR3_VEND','CDR3_NONGERM','CDR3_JEND')
  
  return(split_items)
}

## split cdr3 and add to patient dataframe
all_patients_pre.split = apply(all_patients_pre.df, 1, split_cdr3)
all_patients_split.df = cbind(all_patients_pre.df,
                            data.frame(t(all_patients_pre.split),stringsAsFactors=F))

## save split data
save(all_patients_split.df,file='data/all_patients_split.RData')



