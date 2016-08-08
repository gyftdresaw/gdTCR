
# evaluating amino acid properties of cdr3 partitioned dataset
setwd('~/Documents/Dinner/gd TCR/TCR_analysis 062116/gdTCR/')

load('data/all_patients_split.RData')

### some CDR3 sequence statistics ###
# define amino acid groupings
# -- simply basic,acidic,polar,nonpolar
hydrophobic_AA = c('A','V','L','I','P','F','M','W','G','C')
hydrophilic_AA = c('N','Q','S','T','Y')
positive_AA = c('K','R','H')
negative_AA = c('D','E')

# function for testing number of group representatives
test_AA = function(seq,AA){
  seq_vector = unlist(strsplit(seq,""))
  return(sum(seq_vector %in% AA)/nchar(seq))
}

test_hscales = function(seq,hscales,scale){
  seq_vector = unlist(strsplit(seq,""))
  return(sum(hscales[[scale]][match(seq_vector,hscales$Residue)])/nchar(seq))
}

## more refined hydrophobicity scales
hscales = read.table('data/hydrophobicity.txt',sep='\t',header=T)

# add cdr3 properties to complete data set
get_aa_properties = function(all_patients.df,expand=T){
  cdr3_properties = all_patients.df
  cdr3_properties$CDR3.length = nchar(as.character(all_patients.df$CDR3))
  cdr3_properties$hydrophobic = unlist(lapply(all_patients.df$CDR3_NONGERM,function(x){return(test_AA(x,hydrophobic_AA))}))
  cdr3_properties$hydrophilic = unlist(lapply(all_patients.df$CDR3_NONGERM,function(x){return(test_AA(x,hydrophilic_AA))}))
  cdr3_properties$positive = unlist(lapply(all_patients.df$CDR3_NONGERM,function(x){return(test_AA(x,positive_AA))}))
  cdr3_properties$negative = unlist(lapply(all_patients.df$CDR3_NONGERM,function(x){return(test_AA(x,negative_AA))}))
  ## better hydrophobic scales
  cdr3_properties$KD = unlist(lapply(all_patients.df$CDR3_NONGERM,function(x){return(test_hscales(x,hscales,'KD'))}))
  cdr3_properties$WW = unlist(lapply(all_patients.df$CDR3_NONGERM,function(x){return(test_hscales(x,hscales,'WW'))}))
  cdr3_properties$HH = unlist(lapply(all_patients.df$CDR3_NONGERM,function(x){return(test_hscales(x,hscales,'HH'))}))
  
  # expand by count data so boxplot statistics are correct
  cdr3_prop_expand = cdr3_properties[rep(1:nrow(cdr3_properties),cdr3_properties[['Count']]),]

  if (expand)
    return(cdr3_prop_expand)
  else
    return(cdr3_properties)
}

cdr3_prop = get_aa_properties(all_patients_split.df,expand=F)
cdr3_prop_expand = get_aa_properties(all_patients_split.df,expand=T)

require(ggbeeswarm)

## properties not separated by individual, aggregated by group
plot_aggregate_properties = function(prop.df,properties,flabel) {
  for (i in 1:length(plot_properties)) {
    g = ggplot(prop.df,aes_string('Group',properties[i]))
    g = g + geom_quasirandom(aes(color=Group))
    g = g + geom_boxplot(outlier.size=0,alpha=0,width=0.5)
    g = g + facet_grid(Chain ~ Tissue)
    g
    ggsave(paste('vquest_plots/',properties[i],'_',flabel,'.png',sep=''))
  }
}

# plot all properties
plot_properties = colnames(cdr3_prop)[21:28]
plot_aggregate_properties(cdr3_prop,plot_properties,'unique')
plot_aggregate_properties(cdr3_prop_expand,plot_properties,'expand')

### old style plots -- too cluttered ###
## plot properties by group
plot_folder = 'vquest_plots/'
tissues = c('IEL','PBL')
for (i in 1:length(tissues)){
  # length
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),CDR3.length))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + ylab('CDR3 Length')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'cdr3length_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # hydrophobic
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),hydrophobic))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # hydrophilic
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),hydrophilic))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophilic_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # positive
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),positive))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'positive_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # negative
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),negative))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'negative_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # more refined hydrophobic scales
  # Kyte Doolittle (KD)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),KD))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_KD_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # Whimley White (WW)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),WW))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_WW_',tissues[i],'.png',sep=''),width=10,height=5)
  
  # Hessa von Heigne (HH)
  g = ggplot(cdr3_prop_expand[cdr3_prop_expand$Tissue==tissues[i],],aes(factor(ID),HH))
  g = g + geom_boxplot(aes(color=Chain))
  g = g + facet_grid(~Group,scales='free_x',space='free_x')
  g = g + xlab('Patient ID')
  g = g + scale_color_manual(values=c('TRD'='deeppink','TRG'='slateblue'))
  g = g + theme(text = element_text(size=25))
  g
  ggsave(paste(plot_folder,'hydrophobic_HH_',tissues[i],'.png',sep=''),width=10,height=5)
  
}




