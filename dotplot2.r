#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "res"
}

data <- read.csv(args[1], sep='\t', header=TRUE)
data
data$GI_Score<- as.character(data$GI_Score)
data$Project<- as.character(data$Project)
#Then turn it back into a factor with the levels in the correct order
data$GI_Score<- factor(data$GI_Score, levels=unique(data$GI_Score))
data$Project<- factor(data$Project, levels=unique(data$Project))

library(ggplot2)
dot_plot <-ggplot(data, gridColor="white")
dot_plot + geom_point(aes(x=GI_Score, y=Project,size=log10_corrected_pvalues, color=Difference)) + scale_color_manual(values=c('red','blue','grey')) + theme_linedraw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust= 0.5))

fn1 = strsplit(args[1],'/')
print(fn1)
fn2 = paste(args[2], '/', fn1[[1]], '_dotplot.png', sep = '')
print(fn2)
ggsave(filename = fn2[2],width = 15,height =4,dpi = 450,units = "in",device = 'png')
