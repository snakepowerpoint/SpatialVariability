wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "output\\smap"))

data = read.csv("Smap_abundance_coefficients_vs_b.csv", header = TRUE)
data = na.omit(data)
data$pvalue =  sub(pattern = "<", replacement = "", x = data$pvalue)
data$pvalue = as.numeric(data$pvalue)

data = subset(data, subset = data$pvalue < 0.1)

par(pty='s')
plot(abundance_coeff ~ b, data = data, pch = 19, 
     xlim = c(min(data$b)-0.2, max(data$b)+0.2),
     ylim = c(min(data$abundance_coeff)-0.05, max(data$abundance_coeff)+0.05),
     yaxs='i', xaxs='i',
     xlab = expression(italic('b')), ylab = 'Causal effect of abundance')
abline(h = 0)
abline(v = 2)

col = 127/255
rect(xleft = 2, ybottom = 0, 
     xright = max(data$b)+0.2, ytop = max(data$b)+0.05,
     col= rgb(col,col,col,alpha=0.1), border = 'transparent')
rect(xleft = min(data$b)-0.2, ybottom = min(data$abundance_coeff)-0.05, 
     xright = 2, ytop = 0,
     col= rgb(col,col,col,alpha=0.1), border = 'transparent')

# ggplot
library(ggplot2)

