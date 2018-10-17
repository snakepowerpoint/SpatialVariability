wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
#load(file=paste0(wd, "data\\compiledDataForEDM.RData"))
setwd(paste0(wd, "output"))

library(Rcpp)
library(devtools)
library(rEDM)

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

### Extract standardized data
extractStandard = function(data){
    subdata = subset(data, select = c(Year, Quarter, CV.CPUE.std, Shannon.age.std, Total.CPUE.std,
                                      AMO.std, MeanBT.std, CVofBT.std))
    colnames(subdata) = c('Year', 'Quarter', 'CV.CPUE', 'Shannon.age', 
                          'Total.CPUE', 'AMO', 'Mean.BT', 'CV.BT')
    
    return(subdata)
}

result.ch = extractStandard(result.ch)
result.gm = extractStandard(result.gm)
result.ma = extractStandard(result.ma)
result.mm = extractStandard(result.mm)
result.pp = extractStandard(result.pp[19:48, ])  # drop discontinuous data
result.pv = extractStandard(result.pv)
result.ss = extractStandard(result.ss)
result.ssp = extractStandard(result.ssp)
result.te = extractStandard(result.te)

### a function returning all local peaks
# (will be applied to get all local optimal embedding dimensions)
detectPeak = function(data){
    n = length(data)
    peaks = c()
    values = c()
    for (idx in 2:(n-1)){
        if ((data[idx-1] < data[idx]) & (data[idx+1] < data[idx])){
            peaks = c(peaks, idx)
            values = c(values, data[idx])}
    }
    output = data.frame(peaks, values)
    
    return(output[order(output$values, decreasing = TRUE), ])
}

### a function using simplex projection to determine embedding dimension for each variable
plot.simplex = function(data){
    # create a list to record embedding dimensions
    dim.list = list()
    
    # determine embedding dimension of each variables
    for (column in colnames(data)){
        output <- simplex(data[[column]], E = 1:10, tau = 1, tp = 1) # embedding dimension 
        plot(output$E, output$rho, type = "l", main = column, lwd = 2,
             xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")
        dim.list[[column]] = detectPeak(output$rho)
        #dim.list[[column]] = output[, c("E", "rho")][order(output$rho, decreasing = TRUE), ]
    }
    
    return(dim.list)
}


### a function using CCM to determine causal variables (including lagged terms)
# Spatial CV ~ Shannon age, Total CPUE, AMO, Mean BT, CV of BT
library(data.table)

plot.ccm = function(data, lags=8, species){
    dim.list = plot.simplex(data)
    
    cols = c('length', 'lib.size', 'E', 'library', 'target', 'tar.lag', 
             'rho', 'sd.rho', 'kendall.tau', 'significance')
    cv.var = data.frame(matrix(0, length(-lags:0), length(cols)))  # ccm results (cv xmap var)
    colnames(cv.var) = cols
    
    output = data.frame(matrix(0, ncol = length(cols), nrow = 0))  # ccm results for all variables
    colnames(output) = colnames(cv.var)
    
    num_samples = 100
    for (idx in 2:length(data)){
        for (lag in -lags:0){
            cv_x <- ccm(data, E = dim.list[[idx]]$peaks[1], tau = 1, tp = lag, lib_sizes = seq(5, dim(data)[1], 3),
                        lib_column = 1, target_column = idx, num_samples = num_samples)
            cv_x_m <- data.frame(ccm_means(cv_x), sd.rho = with(cv_x, tapply(rho, lib_size, sd)))
            
            n = dim(cv_x_m)[1]  # last row: ccm results with maximal library size
            
            cv.var[(lag+9), "length"] = dim(data)[1]
            cv.var[(lag+9), "lib.size"] = cv_x_m$lib_size[n]
            cv.var[(lag+9), "E"] = dim.list[[idx]]$peaks[1]
            cv.var[(lag+9), "library"] = "cv"
            cv.var[(lag+9), "target"] = names(dim.list)[idx]
            cv.var[(lag+9), "tar.lag"] = lag
            cv.var[(lag+9), "rho"] = cv_x_m$rho[n]
            cv.var[(lag+9), "sd.rho"] = cv_x_m$sd.rho[n]
            cv.var[(lag+9), "kendall.tau"] = cor.test(cv_x_m$lib_size, cv_x_m$rho, method = 'kendall')$p.value
            cv.var[(lag+9), "significance"] = t.test(subset(cv_x, subset = lib_size==cv_x_m$lib_size[n])$rho,
                                                     alternative = "greater")$p.value
        } 
        
        # plot ccm results
        # cv xmap variables
        plot(x = -lags:0, y = cv.var$rho, type = 'l', col = 'blue', ylim = c(-0.5,1), xaxt = 'n',
             xlab = expression(paste('Cross map lag (', italic('l'), ')')), main = species, 
             ylab = expression(paste('Correlation coefficient ( ', rho, ' )')))
        axis(1, at = seq(-8, 0, 1))
        segments(-lags:0, cv.var[, 'rho'] - cv.var[, 'sd.rho'],
                 -lags:0, cv.var[, 'rho'] + cv.var[, 'sd.rho'], col = 'blue')
        segments(-lags:0 - 0.2, cv.var[, 'rho'] - cv.var[, 'sd.rho'],
                 -lags:0 + 0.2, cv.var[, 'rho'] - cv.var[, 'sd.rho'], col = 'blue')
        segments(-lags:0 - 0.2, cv.var[, 'rho'] + cv.var[, 'sd.rho'],
                 -lags:0 + 0.2, cv.var[, 'rho'] + cv.var[, 'sd.rho'], col = 'blue')
        abline(h = 0)
        legend(x = -3, y = 0.98, legend = paste0('cv xmap ', names(dim.list)[idx]), text.col = c('blue'))
        
        output = rbind(output, cv.var)
    }
    
    output$target = factor(output$target, levels=unique(output$target))
    output = split(output, output$target)
    
    return(c(dim.list[1], output))  #combine embedding dimension of spatial CV with results
}

ccm.ch = plot.ccm(data = result.ch[, -c(1,2)], species = "Clupea harengus")
ccm.gm = plot.ccm(data = result.gm[, -c(1,2)], species = "Gadus morhua")
ccm.ma = plot.ccm(data = result.ma[, -c(1,2)], species = "Melanogrammus aeglefinus")
ccm.mm = plot.ccm(data = result.mm[, -c(1,2)], species = "Merlangius merlangus")
ccm.pp = plot.ccm(data = result.pp[, -c(1,2)], species = "Pleuronectes platessa")
ccm.pv = plot.ccm(data = result.pv[, -c(1,2)], species = "Pollachius virens")
ccm.ss = plot.ccm(data = result.ss[, -c(1,2)], species = "Scomber scombrus")
ccm.ssp = plot.ccm(data = result.ssp[, -c(1,2)], species = "Sprattus sprattus")
ccm.te = plot.ccm(data = result.te[, -c(1,2)], species = "Trisopterus esmarkii")


# keep lags which pass the testing of kendall's tau and t-test
# and find the one having maximal rho
maxSignificant = function(x){
    data = data.frame(x)
    data = subset(data, subset = data$kendall.tau < 0.05 & data$significance < 0.05)
    
    if (nrow(data) > 0){
        return(data[order(data$rho, decreasing = TRUE), ])
    }
    
    else{
        return(NULL)
    }
}

ccm.ch = c(ccm.ch[1], lapply(ccm.ch[-1], FUN = maxSignificant))
ccm.gm = c(ccm.gm[1], lapply(ccm.gm[-1], FUN = maxSignificant))
ccm.ma = c(ccm.ma[1], lapply(ccm.ma[-1], FUN = maxSignificant))
ccm.mm = c(ccm.mm[1], lapply(ccm.mm[-1], FUN = maxSignificant))
ccm.pp = c(ccm.pp[1], lapply(ccm.pp[-1], FUN = maxSignificant))
ccm.pv = c(ccm.pv[1], lapply(ccm.pv[-1], FUN = maxSignificant))
ccm.ss = c(ccm.ss[1], lapply(ccm.ss[-1], FUN = maxSignificant))
ccm.ssp = c(ccm.ssp[1], lapply(ccm.ssp[-1], FUN = maxSignificant))
ccm.te = c(ccm.te[1], lapply(ccm.te[-1], FUN = maxSignificant))

### need to modify saving function
# save CCM results
saveCCM = function(ccm_list, species){
    significant_lagged_var = ccm_list[-1]
    significant_lagged_var = do.call(rbind, significant_lagged_var)
    rownames(significant_lagged_var) = NULL
    significant_lagged_var = subset(significant_lagged_var, 
                                    select=c(target, tar.lag, rho, sd.rho, kendall.tau))
    
    write.csv(significant_lagged_var, paste0(wd, "output\\ccm\\", "ccm_", species, ".csv"))
    
    return(significant_lagged_var)
}

saveCCM(ccm.ch, species = "Clupea harengus")
saveCCM(ccm.gm, species = "Gadus morhua")
saveCCM(ccm.ma, species = "Melanogrammus aeglefinus")
saveCCM(ccm.mm, species = "Merlangius merlangus")
saveCCM(ccm.pp, species = "Pleuronectes platessa")
saveCCM(ccm.pv, species = "Pollachius virens")
saveCCM(ccm.ss, species = "Scomber scombrus")
saveCCM(ccm.ssp, species = "Sprattus sprattus")
saveCCM(ccm.te, species = "Trisopterus esmarkii")


### plot S-map coefficients results
library(ggplot2)

plot.coeff = function(coefficient_data, lags, species=NULL, mode = "series", rho){
    # coefficient_data contains coefficients for all variables excluding CV.CPUE
    ntime = dim(coefficient_data)[1]
    nvar = dim(coefficient_data)[2]
    coeff.melt = cbind(date = rep(1:ntime, nvar), melt(coefficient_data))
    
    whichvar = match(colnames(coefficient_data), order.var[-1])
    cl = rep(c(2:6)[whichvar], each = ntime) #points colour
    #points notation: c('Shannon.age','Total.CPUE','AMO','Mean.BT','CV.BT')
    #circle(s),circle,diamond(s),triangle(s),triangle: c(21,1,23,24,2)
    sh = rep(c(21,1,23,24,2)[whichvar], each = ntime)
    
    scaleFUN = function(x){sprintf("%.2f", x)}
    if (mode == "series"){ 
        #plot time series
        smaptime = ggplot(data = coeff.melt, aes(x = date, y = value, group = variable))
        print(smaptime + geom_point(size = 2, col = cl, shape = sh, fill = cl) +
                  geom_line(col = cl) +
                  geom_hline(yintercept = 0, lty = 2) +
                  theme(plot.title = element_text(hjust = 0.5, size = 20),
                        axis.title = element_text(size = 18, face = "bold"),
                        axis.text = element_text(size = 18, colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
                        panel.background = element_blank()) + 
                  xlab('Time') +
                  ylab('S-map coefficients') +
                  scale_y_continuous(labels=scaleFUN) +
                  ggtitle(bquote(paste(italic(.(species)), ", ", rho, " = ", .(rho)))))
    }
    
    if (mode == "boxplot"){
        #plot box plot
        smapbox = ggplot(data = coeff.melt, aes(x = variable, y = value, group = variable))
        print(smapbox + geom_boxplot(na.rm = T, col = c(2:6)[whichvar], lwd = 1, width = 0.5*nvar/5) + 
                  geom_hline(yintercept = 0, lty = 2) +
                  theme(axis.title.x = element_blank(), 
                        plot.title = element_text(hjust = 0.5, size = 20),
                        axis.title = element_text(size = 18, face = "bold"),
                        axis.text = element_text(size = 18, colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
                        panel.background = element_blank()) + 
                  ylab('S-map coefficients') +
                  scale_y_continuous(labels=scaleFUN) +
                  scale_x_discrete(labels=paste0(colnames(coefficient_data), "-", lags)) + 
                  ggtitle(bquote(paste(italic(.(species)), ", ", rho, " = ", .(rho)))))
    }
}


### a function to perform S-map
plot.smap = function(ccm_list, data, species, mode = "boxplot"){
    dim.cv = ccm_list$CV.CPUE
    variables = sapply(ccm_list[-1], FUN=function(ccm_each_var){ccm_each_var[1, ]}, simplify = FALSE)
    variables = do.call(rbind, variables)
    variables = variables[order(variables$rho, decreasing = TRUE), c("target", "tar.lag", "rho")]
    variables$target = as.character(variables$target)
    
    cat("\nThe most significant lagged variables determined by CCM: \n")
    print(variables)
    
    total.var = dim(variables)[1]
    dim.cv$is.smap.embed = (dim.cv$peaks - 1) <= total.var
    
    for (idx in 1:length(dim.cv$peaks)){
        if (dim.cv$is.smap.embed[idx]){
            nvar = dim.cv$peaks[idx]
            nsample = dim(data)[1]
            
            smap.data = data.frame(matrix(0, ncol = nvar, nrow = nsample))
            colnames(smap.data) = c("CV.CPUE", variables$target[1:(nvar-1)])
            
            smap.data[, "CV.CPUE"] = data[, "CV.CPUE"]
            for (idx_var in 1:(nvar-1)){
                lags = abs(variables$tar.lag[idx_var])
                vars = variables$target[idx_var]
                smap.data[, vars] = c(rep(NA,lags), data[, vars][1:(nsample-lags)])
            }
            #sort smap.data
            order.var <<- c("CV.CPUE", "Shannon.age", "Total.CPUE", "AMO", "Mean.BT", "CV.BT")
            smap.data = smap.data[, order(match(colnames(smap.data), order.var))]
            
            cat("\nData used for S-map: \n")
            print(head(smap.data))
            
            # determine the optimal theta
            theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 
                      0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 
                      2, 3, 4, 6, 8)
            rho.theta = data.frame(matrix(0, length(theta), 2))
            rho.theta[, 1] = theta
            colnames(rho.theta) = c('theta', 'rho')
            
            # compute rho under each theta
            for (i in 1:length(theta)){
                block_lnlp_output <- block_lnlp(smap.data, columns = c(1:dim(smap.data)[2]), target_column = 1, tp = 1,
                                                method = "s-map", num_neighbors = 0, stats_only = F, theta = rho.theta[i, 1],
                                                save_smap_coefficients = T)
                rho.theta[i, "rho"] = block_lnlp_output$rho
            }  # 18 warnings are OK
            theta.opt = rho.theta$theta[which.max(rho.theta$rho)]
            
            block_lnlp_output <- block_lnlp(smap.data, columns = c(1:dim(smap.data)[2]), target_column = 1, tp = 1, 
                                            method = "s-map", num_neighbors = 0, stats_only = F, theta = theta.opt,
                                            save_smap_coefficients = T)
            
            # The first E columns are for the E lags or causal variables,
            # while the (E+1)th column is the constant
            # we plot only the causal variables
            coeff = data.frame(block_lnlp_output$smap_coefficients[[1]])  # s-map coefficient
            colnames(coeff) = c(colnames(smap.data), "Constant")
            
            cat("\nS-map coefficient: \n")
            print(head(coeff))
            
            plot.coeff(coeff[, -c(1, length(coeff)), drop = FALSE], lags = colSums(is.na(smap.data[,-1,drop=FALSE])), 
                       species = species, mode = mode, rho = round(block_lnlp_output$rho, 2))
            ggsave(paste0(wd, "output\\smap\\smap_boxplot_", species, idx, ".png"))
            print(apply(coeff, 2, mean, na.rm = TRUE))
        }
    }
}

plot.smap(ccm.ch, data = result.ch[, -c(1,2)], species = "Clupea harengus")
plot.smap(ccm.gm, data = result.gm[, -c(1,2)], species = "Gadus morhua")
plot.smap(ccm.ma, data = result.ma[, -c(1,2)], species = "Melanogrammus aeglefinus")
plot.smap(ccm.mm, data = result.mm[, -c(1,2)], species = "Merlangius merlangus")
plot.smap(ccm.pp, data = result.pp[, -c(1,2)], species = "Pleuronectes platessa")
plot.smap(ccm.pv, data = result.pv[, -c(1,2)], species = "Pollachius virens")
plot.smap(ccm.ss, data = result.ss[, -c(1,2)], species = "Scomber scombrus")
plot.smap(ccm.ssp, data = result.ssp[, -c(1,2)], species = "Sprattus sprattus")
plot.smap(ccm.te, data = result.te[, -c(1,2)], species = "Trisopterus esmarkii")



# to do list





### Appendix
#species name
"
Clupea harengus
Gadus morhua
Melanogrammus aeglefinus
Merlangius merlangus
Pleuronectes platessa
Pollachius virens
Scomber scombrus
Sprattus sprattus
Trisopterus esmarkii
"
