# save CCM results
saveCCM = function(ccm_list, species){
    rho_for_lags = sapply(ccm_list[-1], FUN = function(ccm_result){return(ccm_result$rho)})
    rho_for_lags = t(rho_for_lags)
    colnames(rho_for_lags) = 8:0
    
    want.var.order = c('Shannon.age','Total.CPUE','AMO','Mean.BT','CV.BT')
    target.var.order = rownames(rho_for_lags)
    rho_for_lags = rho_for_lags[match(want.var.order, target.var.order), ]
    
    write.csv(rho_for_lags, paste0(wd, "output\\", "ccm_", species, ".csv"))
    
    return(rho_for_lags)
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
