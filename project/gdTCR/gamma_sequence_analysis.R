
## More in-depth sequence based analysis of Vgamma chains
## CDR3 non-germline encoded sequences parsed

setwd('~/Documents/Dinner/gd TCR/TCR_analysis-master 062116/gdTCR/')

load('all_patients_gamma.RData')

require(Biostrings)
require(stringr)

## simple amino acid usage ## 
# total counts
aa_counts = t(sapply(all_patients_gamma.df$CDR3_NONGERM,function(x) str_count(x,AA_STANDARD)))
aa_counts.df = data.frame(aa_counts)
colnames(aa_counts.df) = AA_STANDARD
# binary presence/absence
aa_bin.df = aa_counts.df > 0

# combine aa counts with patient data frame
all_aa_counts.df = cbind(all_patients_gamma.df,aa_counts.df)
all_aa_bin.df = cbind(all_patients_gamma.df,aa_bin.df)
# expanded
all_aa_counts_expand.df = all_aa_counts.df[rep(1:nrow(all_aa_counts.df),all_aa_counts.df$Count),]
all_aa_bin_expand.df = all_aa_bin.df[rep(1:nrow(all_aa_bin.df),all_aa_bin.df$Count),]

# histogram plot
require(reshape2)
require(ggplot2)
require(RColorBrewer)

## input is patient dataframe with amino acid counts
plot_aa_usage = function(aa.df,flabel) {
  # aggregate by group and tissue
  aa_summary.df = ddply(aa.df,.(Group,Tissue),colwise(function(x) sum(x)/length(x),AA_STANDARD))
  aa_summary.mlt = melt(aa_summary.df,value.name='freq',variable.name='AA',id.vars=c('Group','Tissue'))
  
  # plot without individual contributions
  g = ggplot(aa_summary.mlt,aes(AA,freq))
  g = g + geom_bar(aes(fill=Group),position='dodge',stat='identity')
  g = g + facet_grid(Tissue ~ .)
  g
  ggsave(paste('vquest_plots/CDR3_nongerm_AA_usage_',flabel,'.png',sep=''))
  
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
  ggsave(paste('vquest_plots/CDR3_nongerm_AA_usage_',flabel,'_indiv.png',sep=''))

}

## input is patient dataframe with amino acid counts
plot_aa_usage_tissue = function(aa.df,flabel) {
  # aggregate by group and tissue
  aa_summary.df = ddply(aa.df,.(Group,Tissue),colwise(function(x) sum(x)/length(x),AA_STANDARD))
  aa_summary.mlt = melt(aa_summary.df,value.name='freq',variable.name='AA',id.vars=c('Group','Tissue'))
  
  # plot without individual contributions
  g = ggplot(aa_summary.mlt,aes(AA,freq))
  g = g + geom_bar(aes(fill=Tissue),position='dodge',stat='identity')
  g = g + facet_grid(Group ~ .)
  g
  ggsave(paste('vquest_plots/CDR3_nongerm_AA_usage_tissue_',flabel,'.png',sep=''))
  
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
  g = ggplot(aa_indiv_summary.mlt,aes(Tissue,freq))
  g = g + geom_bar(aes(fill=factor(ID)),stat='identity')
  g = g + scale_fill_manual(values=indivcolors)
  g = g + facet_grid(Group ~ AA)
  g
  ggsave(paste('vquest_plots/CDR3_nongerm_AA_usage_tissue_',flabel,'_indiv.png',sep=''))
  
}


# plot all types of counting
plot_aa_usage(all_aa_counts.df,'counts')
plot_aa_usage(all_aa_bin.df,'binary')
plot_aa_usage(all_aa_counts_expand.df,'counts_expand')
plot_aa_usage(all_aa_bin_expand.df,'binary_expand')
# plot all types of counting -- tissue comparison
plot_aa_usage_tissue(all_aa_counts.df,'counts')
plot_aa_usage_tissue(all_aa_bin.df,'binary')
plot_aa_usage_tissue(all_aa_counts_expand.df,'counts_expand')
plot_aa_usage_tissue(all_aa_bin_expand.df,'binary_expand')






