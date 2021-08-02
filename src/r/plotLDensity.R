#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)
library(ggplot2)
library(data.table)

data <- read.table('%3%/%4%', sep="\t", header=T)
plotLadderDensity(data$READ, data$COPY, xl="Read Count", yl="Density", title="")

<<@@@@>>
