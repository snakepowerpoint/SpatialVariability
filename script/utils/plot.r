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
    var_order = c('AgeDiversity', 'Abundance')
    data$library = factor(data$library, levels=var_order)
    species = species_list$species
    
    ggplot(aes(x = tar.lag, y = rho, color=library), data = data) + theme_bw() +
        geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black') +
        geom_line(aes(color = library), size = 1) +
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

