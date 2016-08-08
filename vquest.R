
# parsing full vquest output -- excel file with multiple sheets
setwd('~/Documents/Dinner/gd TCR/TCR_analysis-master 062116/gdTCR/')

# require(xlsx)
require(XLConnect) # this package lets you import all sheets at once
require(Biostrings) # for nt translation

# helper function to extract matches from regex
get_matches = function(mstring,mresult){
  matches = rep(NA,length(mresult))
  for (i in 1:length(mresult)){
    matches[i] = substr(mstring,mresult[i],mresult[i]+attributes(mresult)$match.length[i]-1)
  }
  return(matches)
}

# helper to parse recorded genes and alleles
get_genes = function(gene_strings){
  match_results = gregexpr("TR\\S*",gene_strings)
  
  genes.df = data.frame(Gene=character(0),Allele=character(0),stringsAsFactors=F)
  alts = list()
  for (i in 1:length(gene_strings)){
    matches = get_matches(gene_strings[i],match_results[[i]])
    matches.lst = strsplit(matches,'\\*')
    genes.df[i,] = matches.lst[[1]]
    if (length(matches) > 1)
      alts[[i]] = matches.lst[-1]
    else
      alts[[i]] = NA
  }
  
  return(list(genes.df,alts))
}

## try reading one of the vquest spreadsheets
parse_vquest_gamma = function(vfile) {
  # vfile = 'data/vquest/KL304_vquest 2'
  wb = loadWorkbook(vfile)
  tst.lst = readWorksheet(wb,getSheets(wb))
  
  # this is a list of data frames named by the sheet
  summary.df = tst.lst[[which(names(tst.lst)=='Summary')]]
  nt.df = tst.lst[[which(names(tst.lst)=='Junction')]]
  # aa.df = tst.lst[[which(names(tst.lst)=='AA-sequences')]][prod_filter,]
  
  ## limit to number of entries in junction sheet
  summary.df = summary.df[1:nrow(nt.df),]
  
  # unique(summary.df$Functionality)
  # always 'productive','No results', 'unproductive (see comment)', or 'unknown (see comment)'
  # only going to retain productive rearrangements
  prod_filter = summary.df$Functionality == 'productive'

  ## extracting V and J gene
  v_genes = paste(summary.df[prod_filter,]$V.GENE.and.allele,
                  summary.df[prod_filter,]$V.REGION.potential.ins.del,sep=' ')
  j_genes = paste(summary.df[prod_filter,]$J.GENE.and.allele,
                  summary.df[prod_filter,]$J.GENE.and.allele.comment,sep=' ')
  
  # data frames with genes and alleles for each productive sample
  v_genes.df = get_genes(v_genes)[[1]]
  colnames(v_genes.df) = c('VGene','VAllele')
  v_alts = get_genes(v_genes)[[2]]
  j_genes.df = get_genes(j_genes)[[1]]
  colnames(j_genes.df) = c('JGene','JAllele')
  j_alts = get_genes(j_genes)[[2]]
  
  ### individual nucleotide sequences component separated ###
  nt.df = nt.df[prod_filter,]
  
  ## extracting nucleotide sequence components
  # the junction is already separated - don't need to perform end matching
  
  nt_concise.df = nt.df[,c('Sequence.number','Sequence.ID','JUNCTION',
                           'X3.V.REGION','P3.V','N.REGION','P5.J','X5.J.REGION')]
  # replace NA values with empty strings
  nt_concise.df[is.na(nt_concise.df)] = ''
  
  ## need to reshuffle these nt to be codon units
  ## only retaining whole codons in conserved V and J ends
  
  # V end without trailing nucleotides
  nt_concise.df$VEND = sapply(nt_concise.df$X3.V.REGION,
                              function(x) subseq(x,end=-((nchar(x) %% 3)+1)))
  # trailing nucleotides from V end
  nt_concise.df$VNT = sapply(nt_concise.df$X3.V.REGION,
                             function(x) subseq(x,width=nchar(x) %% 3,end=-1))
  # V check
  if (!all(paste(nt_concise.df$VEND,nt_concise.df$VNT,sep='') == nt_concise.df$X3.V.REGION))
    print('V ends not correctly parsed')
  # trailing nucleotides from J end
  nt_concise.df$JNT = sapply(nt_concise.df$X5.J.REGION,
                             function(x) subseq(x,width=nchar(x) %% 3,start=1))
  # J end without trailing nucleotides
  nt_concise.df$JEND = sapply(nt_concise.df$X5.J.REGION,
                             function(x) subseq(x,start=(nchar(x) %% 3)+1))
  # J check
  if (!all(paste(nt_concise.df$JNT,nt_concise.df$JEND,sep='') == nt_concise.df$X5.J.REGION))
    print('J ends not correctly parsed')
  
  # now we just need to concatenate the middle regions
  nt_concise.df$NONGERM = apply(nt_concise.df[,c('VNT','P3.V','N.REGION','P5.J','JNT')],1,
                                  function(x) paste(x,collapse=''))
  # checks
  if (!all(paste(nt_concise.df$VEND,nt_concise.df$NONGERM,nt_concise.df$JEND,sep='') == 
      nt_concise.df$JUNCTION))
      print('Components do not match junction')
  
  # convert nt sequences to AA sequences
  nt_concise.df$JUNCTION_AA = sapply(nt_concise.df$JUNCTION,
                                     function(x) toString(translate(DNAString(x))))
  nt_concise.df$VEND_AA = sapply(nt_concise.df$VEND,
                                 function(x) toString(translate(DNAString(x))))
  nt_concise.df$NONGERM_AA = sapply(nt_concise.df$NONGERM,
                                 function(x) toString(translate(DNAString(x))))
  nt_concise.df$JEND_AA = sapply(nt_concise.df$JEND,
                                 function(x) toString(translate(DNAString(x))))
  if (!all(paste(nt_concise.df$VEND_AA,nt_concise.df$NONGERM_AA,nt_concise.df$JEND_AA,sep='') == 
           nt_concise.df$JUNCTION_AA))
    print('AA Components do not match AA junction')

  parsed_tcr.df = data.frame(TRV=v_genes.df$VGene,TRV_Allele=v_genes.df$VAllele,
                             TRJ=j_genes.df$JGene,TRJ_Allele=j_genes.df$JAllele,
                             CDR3=nt_concise.df$JUNCTION_AA,
                             CDR3_VEND=nt_concise.df$VEND_AA,
                             CDR3_NONGERM=nt_concise.df$NONGERM_AA,
                             CDR3_JEND=nt_concise.df$JEND_AA)
  
  return(parsed_tcr.df)
}

## example single usage
parsed_example.df = parse_vquest_gamma('data/vquest/Active/Chicago #22/Chicago #22_IEL_Vd1_TRG_1_KL126_Dec2014_vquest_1.xls')

################################
### parsing complete dataset ###
################################

# 3 patient groups
patient_groups = c('Control','Active','GFD')
patient_ids = list(c(7,13,40,53,106,110,144,111),c(35,46,47,51,81,112,143,22),c(4,9,28,33,41,43,113,3))

# given patient group and id, return tcr sequencing data table
## parsed according to vquest nucleotide sequence
## ONLY for gamma chain files 
read_patient_gamma = function(group,id) {
  # tables to retrieve per individual
  tissues = c('IEL','PBL')
  
  # want to add group,id,tissue,chain columns to data table
  # explicitly written to keep types consistent
  ## MODIFIED to include alleles and parsed CDR3, excluding NT
  patient.df = data.frame(Group=character(0),ID=numeric(0),Tissue=character(0),Chain=character(0),
                          TRV=character(0),TRV_Allele=character(0),TRJ=character(0),TRJ_Allele=character(0),
                          CDR3=character(0),CDR3_VEND=character(0),CDR3_NONGERM=character(0),CDR3_JEND=character(0),
                          Freq=numeric(0),Count=numeric(0),stringsAsFactors=FALSE)
  
  patient_path = paste('data/vquest/',group,'/Chicago #',id,sep='')
  patient_files = list.files(patient_path)
  
  for (i in 1:length(tissues)) {
    tissue.df = data.frame(matrix(nrow=0,ncol=11))
    ## this finds all files with corresponding tissue and chain
    fids = grep(paste(tissues[i],'_Vd1_TRG',sep=''),patient_files)
    for (j in seq_along(fids)) {
      curr_file = paste(patient_path,'/',patient_files[fids[j]],sep='')
      print(paste('parsing',curr_file))
      parsed.df = parse_vquest_gamma(curr_file)
      
      # add columns and merge with current tissue data table
      tissue.df = rbind(tissue.df,cbind(data.frame(Group=group,ID=id,Tissue=tissues[i],
                                    Chain='TRG',stringsAsFactors=FALSE),parsed.df))
      
    }
    # only proceed if files were parsed from patient
    if (length(fids) > 0) {
      ## tissue table is in expanded format, summarize for unique sequence counts
      concise_tissue.df = ddply(tissue.df,colnames(tissue.df),summarize,
                                Count=length(ID),Freq=Count/nrow(tissue.df))
      # order by count/frequency
      concise_tissue.df = concise_tissue.df[order(concise_tissue.df$Count,decreasing=T),]
      
      # concat with patient dataframe
      patient.df = rbind(patient.df,concise_tissue.df)
    }
  }
  return(patient.df)
}

# mostly works, there are some file warnings -- should check
## iterate through all groups and ids and read patient tables
all_patients_gamma.df = data.frame(matrix(nrow=0,ncol=14))
for (i in 1:length(patient_groups)) {
  for (j in 1:length(patient_ids[[i]])) {
    print(paste(patient_groups[i],patient_ids[[i]][j]))
    all_patients_gamma.df = rbind(all_patients_gamma.df,read_patient_gamma(patient_groups[i],patient_ids[[i]][j]))
  }
}

# save gamma result, takes a couple min to run
save(all_patients_gamma.df,file='all_patients_gamma.RData')








