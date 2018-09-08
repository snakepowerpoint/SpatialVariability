load(file="C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\scripts\\data\\length.RData")

### Taylor's power law
# compute total CPUE and mean and variance of CPUE
taylor <- lapply(cpue.species.drop0, FUN=function(y){
    data = data.frame(y)
    result = with(data, aggregate(x, by=list(Year=Year, Quarter=Quarter), FUN=function(z){var(z)}))
    result1 = with(data, aggregate(x, by=list(Year=Year, Quarter=Quarter), FUN=function(z){mean(z)}))
    result2 = with(data, aggregate(x, by=list(Year=Year, Quarter=Quarter), FUN=sum))
    cbind(result, result1$x, result2$x)
})

# change the column names and sort the data
taylor <- lapply(taylor, FUN=function(y){
    data = data.frame(y)
    colnames(y) = c("Year", "Quarter", "Var.CPUE", "Mean.CPUE", "Total.CPUE")
    y[with(y, order(Year, Quarter)), ]
})

# now we have variance and mean of CPUE, and total CPUE for each cruise for each species
# remove years with 0 mean or variance because log(0) will cause error
taylor <- lapply(taylor, FUN=function(y){
    data = data.frame(y)
    choose = apply(data, 1, FUN=function(row){all(row!=0)})
    data[choose, ]
})

# compute Taylor's exponent
taylor.coeff = lapply(taylor, FUN=function(y){
    data = data.frame(y)
    data = na.omit(data)
    if (nrow(data) >= 5){ # only if data points are larger than 5
        fit = lm(log(Var.CPUE)~log(Mean.CPUE), data=data) # fit Taylor's equation
        output = c(as.numeric(fit$coefficients), summary(fit)$coefficients[2,4],
                   as.numeric(confint(fit, parm = 2, level = 0.95)[1,]))
    }
    else {
        output = rep(NA, 5)
    }
    output = data.frame(matrix(output, ncol=5))
    colnames(output) = c("a", "b", "p_val", "interval_low", "interval_up")
    
    return(output)
})
