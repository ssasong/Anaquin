#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t', stringsAsFactors=FALSE)

plotLadderCopy(data$NAME, data$COPY * data$STOICHOMETRY, log2(data$READ), xl="Copy-number (cn)", yl="Read Count (log2)")
<<@@@@>>
