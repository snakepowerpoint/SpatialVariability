wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "output\\"))

species = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus",
            "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
            "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")

result.ch = read.csv(paste0(species[1], ".csv"), header = TRUE)
result.gm = read.csv(paste0(species[2], ".csv"), header = TRUE)
result.ma = read.csv(paste0(species[3], ".csv"), header = TRUE)
result.mm = read.csv(paste0(species[4], ".csv"), header = TRUE)
result.pp = read.csv(paste0(species[5], ".csv"), header = TRUE)
result.pv = read.csv(paste0(species[6], ".csv"), header = TRUE)
result.ss = read.csv(paste0(species[7], ".csv"), header = TRUE)
result.ssp = read.csv(paste0(species[8], ".csv"), header = TRUE)
result.te = read.csv(paste0(species[9], ".csv"), header = TRUE)

# standarize
result.ch = cbind(result.ch[, c(1,2)], scale(result.ch[, -c(1,2)]))
result.gm = cbind(result.gm[, c(1,2)], scale(result.gm[, -c(1,2)]))
result.ma = cbind(result.ma[, c(1,2)], scale(result.ma[, -c(1,2)]))
result.mm = cbind(result.mm[, c(1,2)], scale(result.mm[, -c(1,2)]))
result.pp = cbind(result.pp[, c(1,2)], scale(result.pp[, -c(1,2)]))
result.pv = cbind(result.pv[, c(1,2)], scale(result.pv[, -c(1,2)]))
result.ss = cbind(result.ss[, c(1,2)], scale(result.ss[, -c(1,2)]))
result.ssp = cbind(result.ssp[, c(1,2)], scale(result.ssp[, -c(1,2)]))
result.te = cbind(result.te[, c(1,2)], scale(result.ch[, -c(1,2)]))

# plot BT versus SST, and add regression line
library(ggplot2)
library(devtools)

lm_eqn <- function(df){
    m <- lm(MeanBT ~ MeanSST, df);
    coeff = round(m$coefficients, 2)
    eq <- substitute(italic(y) == a + b ~italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = as.numeric(coeff[1]), 
                          b = as.numeric(coeff[2]), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

# whole data
data = result.ch
ggplot(data, mapping = aes(x = MeanSST, y = MeanBT)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y~x) +
    geom_text(aes(x = -Inf, y = Inf, label = lm_eqn(data)), parse = TRUE, 
              hjust = -0.1, vjust = 2, size = 5) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 18, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          panel.background = element_blank()) +
    xlab('Sea surface temperature') +
    ylab('Bottom temperature')
cor.test(data$MeanSST, data$MeanBT)

# quarter 1
data = subset(result.ch, subset = result.ch$Quarter == 1)
ggplot(data, mapping = aes(x = MeanSST, y =MeanBT)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y~x) +
    geom_text(aes(x = -Inf, y = Inf, label = lm_eqn(data)), parse = TRUE, 
              hjust = -0.1, vjust = 2, size = 5) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 18, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          panel.background = element_blank()) +
    xlab('Sea surface temperature') +
    ylab('Sea bottom temperature')
cor.test(data$MeanSST, data$MeanBT)


# quarter 3
data = subset(result.ch, subset = result.ch$Quarter == 3)
ggplot(data, mapping = aes(x = MeanSST, y =MeanBT)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y~x) +
    geom_text(aes(x = -Inf, y = Inf, label = lm_eqn(data)), parse = TRUE, 
              hjust = -0.1, vjust = 2, size = 5) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 18, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          panel.background = element_blank()) +
    xlab('Sea surface temperature') +
    ylab('Sea bottom temperature')
cor.test(data$MeanSST, data$MeanBT)





### Legacy, please ignore the following codes
library(reshape2)
melt.ch = melt(result.ch, id.vars = "Year")
