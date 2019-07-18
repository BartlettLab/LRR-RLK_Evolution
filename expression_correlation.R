library(ggplot2)
setwd("/Users/jarrettman/Downloads")



tbl <- read.csv("expression.csv")
tbl <- data.frame(tbl[0:140,]) #remove text at bottom


#SET COLUMN of genes you want to compare
y_axis_gene <- 2
x_axis_gene <- 3

#convert to numberic vectors
tbl[,x_axis_gene] <- as.numeric(as.character(tbl[,x_axis_gene]))
tbl[,y_axis_gene] <- as.numeric(as.character(tbl[,y_axis_gene]))

#correlation and plot
cor.test(tbl[,x_axis_gene], tbl[,y_axis_gene], method = "pearson")
ggplot(data = tbl, aes(x = tbl[,x_axis_gene], y = tbl[,y_axis_gene])) + geom_point() + stat_smooth(method = "lm", geom = "smooth") 

