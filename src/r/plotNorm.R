#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv("%3%/%4%", sep="\t")
data <- data[data$S1 != 0 & data$S2 != 0,]

x <- %5%
y <- %6%

plotLinear(row.names(data), x, y, title="%7%", xl="%8%", yl="%9%", showLOQ=FALSE, shouldLog=FALSE, showCol=FALSE)<<@@@@>>