
## More in-depth sequence based analysis of both Vgamma AND Vdelta chains
## CDR3 non-germline encoded sequences parsed

setwd('~/Documents/Dinner/gd TCR/TCR_analysis 062116/gdTCR/')

load('data/all_patients_split.RData')

# separate chains -- only keeping active control gfd
# all_patients_gamma.df = all_patients_split.df[all_patients_split.df$Chain=='TRG' & 
#                                                 all_patients_split.df$Group %in% c('Active','Control','GFD'),]
# all_patients_delta.df = all_patients_split.df[all_patients_split.df$Chain=='TRD' & 
#                                                 all_patients_split.df$Group %in% c('Active','Control','GFD'),]
# only keep active, control, gfd
all_patients_split.df = all_patients_split.df[all_patients_split.df$Group %in% c('Active','Control','GFD'),]

require(Biostrings)
require(stringr)

## simple amino acid usage ## 

# helper for retrieving amino acid counts
get_aa_counts = function(sequences,binary=F) {
  # total counts
  aa_counts = t(sapply(sequences,function(x) str_count(x,AA_STANDARD)))
  aa_counts.df = data.frame(aa_counts)
  colnames(aa_counts.df) = AA_STANDARD
  # binary presence/absence
  aa_bin.df = aa_counts.df > 0
  
  if (binary == T)
    return(aa_bin.df)
  else
    return(aa_counts.df)
}

# helper for expanding dataframe by count column
expand_df = function(a.df){
  return(a.df[rep(1:nrow(a.df),a.df$Count),])
}

## all aa counts combined with data frame
all_counts.df = cbind(all_patients_split.df,get_aa_counts(all_patients_split.df$CDR3_NONGERM,binary=F))
all_bin.df = cbind(all_patients_split.df,get_aa_counts(all_patients_split.df$CDR3_NONGERM,binary=T))
# expanded
all_counts_expand.df = expand_df(all_counts.df)
all_bin_expand.df = expand_df(all_bin.df)

# histogram plot
require(reshape2)
require(ggplot2)
require(RColorBrewer)

## input is patient dataframe with amino acid counts
plot_aa_usage = function(aa_full.df,flabel,chains=c('TRG','TRD'),prelabel='') {
  # for each chain
  for (i in 1:length(chains)) {
    aa.df = aa_full.df[aa_full.df$Chain==chains[i],] # get chain
    
    # aggregate by group and tissue
    aa_summary.df = ddply(aa.df,.(Group,Tissue),colwise(function(x) sum(x)/length(x),AA_STANDARD))
    aa_summary.mlt = melt(aa_summary.df,value.name='freq',variable.name='AA',id.vars=c('Group','Tissue'))
    
    # plot without individual contributions
    g = ggplot(aa_summary.mlt,aes(AA,freq))
    g = g + geom_bar(aes(fill=Group),position='dodge',stat='identity')
    g = g + facet_grid(Tissue ~ .)
    g
    ggsave(paste('vquest_plots/',prelabel,'CDR3_nongerm_AA_usage_',chains[i],'_',flabel,'.png',sep=''))
    
    ## TISSUE ##
    # plot without individual contributions
    g = ggplot(aa_summary.mlt,aes(AA,freq))
    g = g + geom_bar(aes(fill=Tissue),position='dodge',stat='identity')
    g = g + facet_grid(Group ~ .)
    g
    ggsave(paste('vquest_plots/',prelabel,'CDR3_nongerm_AA_usage_tissue_',chains[i],'_',flabel,'.png',sep=''))
    
    # individual contributions displayed
    aa_group_norm = function(t.df){
      glength = sum(aa.df$Group == t.df$Group[1] & aa.df$Tissue == t.df$Tissue[1])
      return(cbind(data.frame(ID=t.df$ID),t.df[,AA_STANDARD]/glength))
    }
    aa_indiv_summary.df = ddply(aa.df,.(Group,Tissue,ID),colwise(function(x) sum(x),AA_STANDARD))
    aa_indiv_norm.df = ddply(aa_indiv_summary.df,.(Group,Tissue),aa_group_norm)
    
    aa_indiv_summary.mlt = melt(aa_indiv_norm.df,value.name='freq',variable.name='AA',id.vars=c('Group','Tissue','ID'))
    
    indivcolors = sample(colorRampPalette(brewer.pal(8,'Set2'))(length(unique(aa_indiv_summary.mlt$ID))))
    # plot with individual contributions
    g = ggplot(aa_indiv_summary.mlt,aes(Group,freq))
    g = g + geom_bar(aes(fill=factor(ID)),stat='identity')
    g = g + scale_fill_manual(values=indivcolors)
    g = g + facet_grid(Tissue ~ AA)
    g
    ggsave(paste('vquest_plots/',prelabel,'CDR3_nongerm_AA_usage_',chains[i],'_',flabel,'_indiv.png',sep=''))
    
    ## TISSUE ##
    # plot with individual contributions
    g = ggplot(aa_indiv_summary.mlt,aes(Tissue,freq))
    g = g + geom_bar(aes(fill=factor(ID)),stat='identity')
    g = g + scale_fill_manual(values=indivcolors)
    g = g + facet_grid(Group ~ AA)
    g
    ggsave(paste('vquest_plots/',prelabel,'CDR3_nongerm_AA_usage_tissue_',chains[i],'_',flabel,'_indiv.png',sep=''))
  }
}

# plot all types of counting
plot_aa_usage(all_counts.df,'counts')
plot_aa_usage(all_bin.df,'binary')
plot_aa_usage(all_counts_expand.df,'counts_expand')
plot_aa_usage(all_bin_expand.df,'binary_expand')


### plot aa usage adjacent to J segment ###
CDR3_NONGERM_LAST = sapply(all_patients_split.df$CDR3_NONGERM,function(x) substr(x,nchar(x),nchar(x)))
j_counts.df = cbind(all_patients_split.df,get_aa_counts(CDR3_NONGERM_LAST,binary=F))
j_bin.df = cbind(all_patients_split.df,get_aa_counts(CDR3_NONGERM_LAST,binary=T))
# expanded
j_counts_expand.df = expand_df(all_counts.df)
j_bin_expand.df = expand_df(all_bin.df)

## plotting usage distribution
plot_aa_usage(j_counts.df,'counts',prelabel='Jadj_')
plot_aa_usage(j_bin.df,'binary',prelabel='Jadj_')
plot_aa_usage(j_counts_expand.df,'counts_expand',prelabel='Jadj_')
plot_aa_usage(j_bin_expand.df,'binary_expand',prelabel='Jadj_')



