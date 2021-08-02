#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv("%3%/%4%", sep="\t")
data <- data[data$LABEL == "Somatic",]

x <- data$EXP_FREQ
y <- data$OBS_FREQ

plotLinear(row.names(data), x, y, title="Allele Frequency Ladder", xl="Expected Allele Frequency (log2)", yl="Measured Allele Frequency (log2)", showLOQ=F, shouldLog=T)<<@@@@>>