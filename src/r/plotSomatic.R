#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv("%3%/%4%", sep="\t")
data <- data[%9% > 0 & !is.na(%9%),]
data <- data[data$LABEL == "Somatic",]

x <- %8%
y <- %9%

plotLinear(row.names(data), x, y, title="%5%", xl="%6%", yl="%7%", showLOQ=F, shouldLog=T)<<@@@@>>