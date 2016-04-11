
# parsing full vquest output -- excel file with multiple sheets
setwd('~/Documents/Dinner/gd TCR/')

# require(xlsx)
require(XLConnect) # this package lets you import all sheets at once

## try reading one of the vquest spreadsheets
wb = loadWorkbook('data/vquest/KL304_vquest 2')
tst.lst = readWorksheet(wb,getSheets(wb))

# this is a list of data frames named by the sheet
summary.df = tst.lst[[which(names(tst.lst)=='Summary')]]

# unique(summary.df$Functionality)
# always 'productive','No results', 'unproductive (see comment)', or 'unknown (see comment)'


