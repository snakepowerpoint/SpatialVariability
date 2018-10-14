wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
#load(file=paste0(wd, "data\\cpue.RData")
setwd(paste0(wd, "data\\"))

cpue = read.csv("compiled_cpue_data.csv", header = TRUE)
cpue = split(cpue, f = cpue$Species)

### Taylor's power law
# compute total CPUE and mean and variance of CPUE
taylor <- lapply(cpue, FUN=function(y){
    data = data.frame(y)
    cpue.var = with(data, aggregate(CPUE, by=list(Year=Year, Quarter=Quarter), FUN=function(z){var(z)}))
    cpue.mean = with(data, aggregate(CPUE, by=list(Year=Year, Quarter=Quarter), FUN=function(z){mean(z)}))
    cpue.sum = with(data, aggregate(CPUE, by=list(Year=Year, Quarter=Quarter), FUN=sum))
    
    out = cbind(cpue.var, cpue.mean$x, cpue.sum$x)
    
    # change the column names and sort the data
    colnames(out) = c("Year", "Quarter", "Var.CPUE", "Mean.CPUE", "Total.CPUE")
    out[with(out, order(Year, Quarter)), ]
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

taylor.dataframe = do.call(rbind, taylor.coeff)
taylor.dataframe$Species = rownames(taylor.dataframe)
rownames(taylor.dataframe) = NULL
taylor.dataframe = cbind(Species = taylor.dataframe$Species, 
                         taylor.dataframe[, -dim(taylor.dataframe)[2]])

write.csv(taylor.dataframe, file=paste0(wd, "output\\taylor.csv"), row.names=FALSE)


#save.image(file="C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\data\\taylor.RData")
