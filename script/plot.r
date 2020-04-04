## plot S-map coefficients
plot_mode = "box"  # 'series' or 'box'

setwd(paste0(wd, 'script'))
source("utils/plot.r")

lapply(EDM_lib_var, function(item, colors=cl, shapes=sh, mode=plot_mode){
    if (is.na(item$Smap1)[1]){
        return(NULL)
    }
    
    species = item$species
    smap_model = grep("Smap", names(item))
    num_smap_model = length(smap_model)
    for (i in 1:num_smap_model){
        smapplot = plotSmapCoeff(smap_result_list=item[[smap_model[i]]],
                                 species=species,
                                 colors=colors,
                                 shapes=shapes,
                                 mode=mode)
        save_path = paste0(wd, smap_path, mode, "_", species, i)
        
        file_name_eps = paste0(save_path, ".eps")
        ggsave(filename=file_name_eps, plot=smapplot, width=9, height=6, units="in")
        file_name_png = paste0(save_path, ".png")
        ggsave(filename=file_name_png, plot=smapplot, width=9, height=6, units="in")
    }
})

if (plot_mode == "series"){
    smap_timeseries_legend(lib_var=library_var, colors=cl, shapes=sh)
} else if (plot_mode == "box"){
    smap_boxplot_legend(lib_var=library_var, colors=cl)
}

file_name_eps = paste0(wd, smap_path, plot_mode, "_legend", ".eps")
ggsave(filename = file_name_eps)
