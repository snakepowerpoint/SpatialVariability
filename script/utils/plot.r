library(ggplot2)
library(extrafont)
library(sysfonts)
library(showtext)

loadfonts(device = "win")
# add the Arial font
font_add("Arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


# function to plot CCM results
plot_ccm_result = function(ccm_result_data, variable, species){
    # AgeDiversity/Abundance xmap fishing mortality
    subdata = subset(ccm_result_data, subset = ccm_result_data$library == variable)
    max_lag = max(abs(subdata$tar.lag))
    
    plot(x = -max_lag:0, y = subdata$rho, type = 'l', col = 'blue', ylim = c(-0.5,1), xaxt = 'n',
         xlab = expression(paste('Cross map lag (', italic('l'), ')')), main = species, 
         ylab = expression(paste('Correlation coefficient ( ', rho, ' )')))
    axis(1, at = seq(-max_lag, 0, 1))
    segments(-max_lag:0, subdata[, 'rho'] - subdata[, 'sd.rho'],
             -max_lag:0, subdata[, 'rho'] + subdata[, 'sd.rho'], col = 'blue')
    segments(-max_lag:0 - 0.1, subdata[, 'rho'] - subdata[, 'sd.rho'],
             -max_lag:0 + 0.1, subdata[, 'rho'] - subdata[, 'sd.rho'], col = 'blue')
    segments(-max_lag:0 - 0.1, subdata[, 'rho'] + subdata[, 'sd.rho'],
             -max_lag:0 + 0.1, subdata[, 'rho'] + subdata[, 'sd.rho'], col = 'blue')
    abline(h = 0)
    legend(x = -max_lag, y = 0.98, legend = paste0(variable, ' xmap fishingM'), text.col = c('blue'))
}

# ggplot function to plot CCM results
gplot_ccm_result = function(species_list){
    # AgeDiversity/Abundance xmap fishing mortality
    data = species_list$ccm
    data[data$kendall.tau >= 0.1 | data$significance >= 0.1, c('rho')] = NA
    var_order = c('AgeDiversity', 'Abundance')
    data$library = factor(data$library, levels=var_order)
    species = species_list$species
    
    ggplot(aes(x = tar.lag, y = rho, color=library), data = data) + theme_bw() +
        geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black') +
        geom_point(size = 5) +
        geom_segment(aes(x=tar.lag-0.1, y=rho+sd.rho, xend=tar.lag+0.1, yend=rho+sd.rho), size=1) +
        geom_segment(aes(x=tar.lag-0.1, y=rho-sd.rho, xend=tar.lag+0.1, yend=rho-sd.rho), size=1) +
        geom_segment(aes(x=tar.lag, y=rho-sd.rho, xend=tar.lag, yend=rho+sd.rho), size=1) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(size = 1.1),
              axis.line = element_line(color = 'black'), 
              axis.title = element_text(color = 'black', size = 16),
              axis.text = element_text(color = 'black', size = 14),
              plot.title = element_text(hjust = 0.5, vjust = 3, size = 18, face = 'bold.italic'),
              legend.position = c(0.025, 0.9),
              legend.background = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 14), 
              legend.justification = c(0, 0)) +
        scale_y_continuous(limits = c(-1, 1.2)) +
        scale_color_manual(labels = paste0(var_order, ' xmap F'), values = c('red', 'green')) + 
        labs(x = 'Lag of fishing mortality', y = expression(rho), title = species)
}

# plot S-map coefficients
plotSmapCoeff = function(
    smap_result_list, 
    species,
    colors,
    shapes,
    mode="series")
{
    # extract lags for each variable
    data_for_smap = smap_result_list$data
    library_var = names(data_for_smap)[1]
    
    target_vars = names(data_for_smap)[-1]
    num_na = colSums(is.na(data_for_smap))
    lag_of_var = as.numeric(num_na[target_vars])
    
    # coefficients of S-map model without library variable and constant
    data_of_coeff = smap_result_list$coefficients
    data_of_coeff = data_of_coeff[target_vars] 
    
    rho = round(smap_result_list$rho, 2)
    
    # sort data according to customized order of variables
    order_var = variables[sort(match(names(data_of_coeff), variables))]
    lag_of_var = lag_of_var[order(match(names(data_of_coeff), order_var))]
    data_of_coeff = data_of_coeff[, order_var, drop = FALSE]
    
    ntime = dim(data_of_coeff)[1]
    nvar = dim(data_of_coeff)[2]
    coeff.melt = cbind(date = rep(1:ntime, nvar), melt(data_of_coeff))
    
    whichvar = match(names(data_of_coeff), variables)
    cl = rep(colors[whichvar], each = ntime) 
    sh = rep(shapes[whichvar], each = ntime)
    
    if (mode == "series"){
        return(smap_timeseries_plot(data=coeff.melt, 
                                    lib_var=library_var,
                                    cl=cl, 
                                    sh=sh, 
                                    species=species, 
                                    rho=rho))
    } else if (mode == "box"){
        return(smap_boxplot(data=coeff.melt, 
                            lib_var=library_var, 
                            cl=cl, 
                            sh=sh,
                            n_var=length(unique(colors)), 
                            species=species, 
                            rho=rho))
    } else {
        stop("mode must be either 'series' or 'box'")
    }
}

# plot time series
smap_timeseries_plot = function(data, lib_var, cl, sh, species, rho){
    scaleFUN = function(x){sprintf("%.2f", x)}
    
    smaptime = 
        ggplot(data=data, aes(x=date, y=value, shape=variable, color=variable, fill=variable)) + 
        geom_point(size=4) +
        geom_line(size=1) +
        geom_hline(yintercept=0, linetype='dashed') +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 28),
              axis.title = element_text(size = 24, face = "bold"),
              axis.text = element_text(size = 20, colour = "black"),
              panel.border = element_rect(size = 1.1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none") + 
        labs(x = 'Time', y = 'S-map coefficients') +
        ggtitle(bquote(paste(italic(.(species)), ", ", rho, " = ", .(rho)))) +
        scale_y_continuous(labels=scaleFUN) +
        scale_shape_manual(values=unique(sh)) +
        scale_color_manual(values=unique(cl)) +
        scale_fill_manual(values=unique(cl))
        
    return(smaptime)
}

# legend for time series plot
smap_timeseries_legend = function(lib_var, colors, shapes){
    data = data.frame(cbind(colors, shapes))
    data$x = 0
    data$y = 0
    data$variable = factor(c(1:dim(data)[1]))
    
    legend_labels = paste0(variables, " effect on ", lib_var)
    
    smaptime = 
        ggplot(data=data, aes(x=x, y=y, shape=variable, color=variable, fill=variable)) + 
        geom_point(size=4) +
        geom_line() +
        theme_bw() + 
        theme(legend.background = element_blank(),
              legend.text = element_text(size = 16),
              legend.key.size = unit(1, 'cm')) + 
        scale_shape_manual(labels=legend_labels, values=shapes) +
        scale_color_manual(labels=legend_labels, values=colors) +
        scale_fill_manual(labels=legend_labels, values=colors)
    
    return(smaptime)
}

# plot box plot
smap_boxplot = function(data, lib_var, cl, sh, n_var, species, rho){
    n = length(unique(data$variable))
    scaleFUN = function(x){sprintf("%.2f", x)}
    
    smapbox = 
        ggplot(data=data, aes(x=variable, y=value, color=variable)) + 
        geom_boxplot(na.rm=T, lwd=1, width=0.7*n/n_var) + 
        geom_hline(yintercept=0, linetype='dashed') +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              axis.title = element_text(size = 18, face = "bold"),
              axis.title.x = element_blank(), 
              axis.text = element_text(size = 18, colour = "black"),
              panel.border = element_rect(size = 1.1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none") + 
        labs(y = 'S-map coefficients') +
        ggtitle(bquote(paste(italic(.(species)), ", ", rho, " = ", .(rho)))) + 
        scale_y_continuous(labels=scaleFUN) +
        scale_x_discrete(labels=unique(data$variable)) + 
        scale_color_manual(values=unique(cl))
    
    return(smapbox)
}

# legend for box plot
smap_boxplot_legend = function(lib_var, colors){
    data = data.frame(colors)
    data$y = 0
    data$variable = factor(c(1:dim(data)[1]))
    
    legend_labels = paste0(variables, " effect on ", lib_var)
    
    smapbox = 
        ggplot(data=data, aes(x=variable, y=y, color=variable)) + 
        geom_boxplot(na.rm=T, lwd=1) + 
        theme_bw() + 
        theme(legend.background = element_blank(),
              legend.text = element_text(size = 16),
              legend.key.size = unit(1, 'cm')) + 
        scale_color_manual(labels=legend_labels, values=colors)
    
    return(smapbox)
}

