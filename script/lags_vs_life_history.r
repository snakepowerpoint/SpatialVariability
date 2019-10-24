wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
lags = read.csv(file = paste0(wd, "output\\ccm\\", "causalVarLag.csv"), header = TRUE)

K = c(NA, 0.23, 0.19, 0.32, 0.11, 0.07, NA, NA, 0.52)
A50 = c(3, 3.8, 2.5, 1.5, 2.5, 4.6, 1.5, 1.5, 2.3)
#A50 = c(2.5, 3.79, 3.96, 1.46, 2.5, 4.6, 1.5, 1.5, 2.3)

lags$K = rep(K, times = as.numeric(table(lags$species)))
lags$A50 = rep(A50, times = as.numeric(table(lags$species)))
lags
lags = subset(lags, subset = !(lags$species %in% c("Clupea harengus", "Scomber scombrus", "Sprattus sprattus")))

## plot function
library(ggplot2)

scatterplot = function(x, y, xlab, ylab, main, quadratic = FALSE){
    data = data.frame(x = x, y = y)
    if (quadratic){
        fit = lm(y ~ x + I(x^2))
    }
    else {
        fit = lm(y ~ x)
    }
    fit.sum = summary(fit)
    fscore = as.numeric(fit.sum$fstatistic)
    fvalue = round(1 - pf(fscore[1], fscore[2], fscore[3]), 4)
    line.type = ifelse(fvalue < 0.1, 1, 2)
    
    prd = data.frame(x = seq(min(x), max(x), length.out = 50))
    err <- predict(fit, newdata = prd, se.fit = TRUE)
    
    prd$lci <- err$fit - 1.96 * err$se.fit
    prd$fit <- err$fit
    prd$uci <- err$fit + 1.96 * err$se.fit
    
    ggplot(prd, aes(x = x, y = fit)) +
        geom_line(linetype = line.type) +
        geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
        geom_point(data = data, aes(x = x, y = y)) +
        ggtitle(bquote(paste(.(main), ", ", italic(p), "-", value, "=", .(fvalue)))) + 
        xlab(xlab) +
        ylab(ylab) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text = element_text(size = 18, colour = "black"),
              panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
              panel.background = element_blank())
}


## A50
data = lags
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "All causal variables"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "All causal variables", "TRUE"))

#biological
data = subset(lags, subset = lags$target %in% c("AgeDiversity", "Abundance"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "Biotic variables"))

#age diversity
data = subset(lags, subset = lags$target %in% c("AgeDiversity"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "Age diversity"))

#Abundance
data = subset(lags, subset = lags$target %in% c("Abundance"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "Abundance"))

#abiotic
data = subset(lags, subset = lags$target %in% c("AMO", "SBT", "SST", "CVofSBT", "CVofSST"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "Abiotic variables"))

#temperature
data = subset(lags, subset = lags$target %in% c("SBT", "SST"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "Temperature"))

#AMO
data = subset(lags, subset = lags$target %in% c("AMO"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "AMO"))

#CV of temperature
data = subset(lags, subset = lags$target %in% c("CVofSBT", "CVofSST"))
with(data, scatterplot(A50, abs(tar.lag)/2, "A50", "Lag", "CV of Temperature"))

#mean each species
john = with(lags, aggregate(tar.lag, by = list(species = species, A50 = A50), FUN = mean))
plot(abs(x)/2 ~ A50, data = john, pch = 19)
fit = lm(abs(x)/2 ~ A50, data = john)
summary(fit)
abline(fit, lty = 2)

with(john, scatterplot(A50, abs(x)/2, "A50", "Lag", "Mean of all variables"))


## K
data = lags
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "All causal variables", quadratic = TRUE))

#biological
data = subset(lags, subset = lags$target %in% c("AgeDiversity", "Abundance"))
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "Biotic variables", quadratic = TRUE))

#age diversity
data = subset(lags, subset = lags$target %in% c("AgeDiversity"))
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "Age diversity", quadratic = TRUE))

#Abundance
data = subset(lags, subset = lags$target %in% c("Abundance"))
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "Abundance", quadratic = TRUE))

#abiotic
data = subset(lags, subset = lags$target %in% c("AMO", "SBT", "SST", "CVofSBT", "CVofSST"))
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "Abiotic variables", quadratic = TRUE))

#temperature
data = subset(lags, subset = lags$target %in% c("SBT", "SST"))
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "Temperature", quadratic = TRUE))

#amo
data = subset(lags, subset = lags$target %in% c("AMO"))
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "AMO", quadratic = TRUE))

#cv of temperature
data = subset(lags, subset = lags$target %in% c("CVofSBT", "CVofSST"))
with(data, scatterplot(K, abs(tar.lag)/2, "K", "Lag", "CV of Temperature", quadratic = TRUE))

#mean each species
data = with(lags, aggregate(tar.lag, by = list(species = species, K = K), FUN = mean))
plot(abs(x)/2 ~ K, data = data, pch = 19)
fit = lm(abs(x)/2 ~ K + I(K^2), data = data)
summary(fit)

with(data, scatterplot(K, abs(x)/2, "K", "Lag", "Mean of all variables", quadratic = TRUE))



##### linear regression
#biological
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity", "Abundance")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity", "Abundance")))
summary(fit)

#age diversity
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity")))
summary(fit)

#Abundance
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("Abundance")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("Abundance")))
summary(fit)

#abiotic
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AMO", "SBT", "SST", "CVofSBT", "CVofSST")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AMO", "SBT", "SST", "CVofSBT", "CVofSST")))
summary(fit)

#temperature
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("SBT", "SST")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("SBT", "SST")))
summary(fit)

#mean each species
data = with(lags, aggregate(tar.lag, by = list(species = species, K = K), FUN = mean))
plot(abs(x)/2 ~ K, data = data)
fit = lm(abs(x)/2 ~ K, data = data)
summary(fit)


## K, polynomial 
plot(abs(tar.lag)/2 ~ K, data = lags)
fit = lm(abs(tar.lag)/2 ~ K + I(K^2), data = lags)
summary(fit)

#biological
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity", "Abundance")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity", "Abundance")))
summary(fit)

#age diversity
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AgeDiversity")))
summary(fit)

#Abundance
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("Abundance")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("Abundance")))
summary(fit)

#abiotic
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AMO", "SBT", "SST", "CVofSBT", "CVofSST")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("AMO", "SBT", "SST", "CVofSBT", "CVofSST")))
summary(fit)

#temperature
plot(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("SBT", "SST")))
fit = lm(abs(tar.lag)/2 ~ K, data = subset(lags, subset = lags$target %in% c("SBT", "SST")))
summary(fit)

#mean each species
data = with(lags, aggregate(tar.lag, by = list(species = species, K = K), FUN = mean))
plot(abs(x)/2 ~ K, data = data)
fit = lm(abs(x)/2 ~ K, data = data)
summary(fit)



##### Appendix
lags = lapply(EDM_lib_var, function(x){
    data = x$ccm_most_significant
    data$species = x$species
    rownames(data) = NULL
    
    return(data)
})

lags = do.call(rbind, lags)
rownames(lags) = NULL  #remove rownames
lags = lags[, c("species", "target", "tar.lag", "count", "rho")]  #sort
lags

write.csv(lags, file = paste0(wd, "output\\ccm\\", "causalVarLag.csv"), row.names = FALSE)

