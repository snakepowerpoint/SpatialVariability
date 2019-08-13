library(gtools)

# permutation function
permu <- function(n, m){
    t(combn(1:n, m, FUN=function(y){
        replace(rep(0,n), y, 1)
    }))
}
#permu(4,2)

# robustness test function
robust = function(data, dim_lib_var, lags){
    n = NROW(data) # length of time series
    num_var = NCOL(data[, -1]) # number of variables, first one is library variable
    lib_var = names(data)[1]
    
    # a mask data frame to choose variables
    mask = data.frame(permu(num_var, dim_lib_var-1)) 
    coeff.mean = apply(mask, 1, FUN=function(x){
        vars = names(data)[-1][x==1] # variables chosen in this permutation
        lags = lags[x==1] # lags of variables
        
        # generate lagged variable
        block.m = subset(data, select = vars)
        block.m = apply(matrix(seq_along(lags)), 1, FUN=function(y){
            if (lags[y] > 0){
                c(rep(NA,lags[y]), block.m[, y][-c((n-lags[y]+1):n)])
            } else {block.m[, y]}
        }) 
        block.m = cbind(subset(data, select=lib_var), block.m)
        
        # find the optimal theta
        theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 
                  0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
        rho.theta = data.frame(matrix(0, length(theta), 2))
        colnames(rho.theta) = c('theta', 'rho')
        rho.theta[, 'theta'] = theta
        
        # compute rho under each theta
        for (i in 1:length(theta)){
            block_lnlp_output <- block_lnlp(block.m, columns = c(1:dim(block.m)[2]), target_column = 1, 
                                            tp = 1, method = "s-map", num_neighbors = 0, stats_only = F, 
                                            theta = rho.theta[i, 'theta'], save_smap_coefficients = T)
            rho.theta[i, 'rho'] = block_lnlp_output$rho
        }
        theta.opt = rho.theta$theta[which.max(rho.theta$rho)]
        
        block_lnlp_output <- block_lnlp(block.m, columns = c(1:dim(block.m)[2]), target_column = 1, 
                                        tp = 1, method = "s-map", num_neighbors = 0, stats_only = F, 
                                        theta = theta.opt, save_smap_coefficients = T)
        coeff = data.frame(block_lnlp_output$smap_coefficients[[1]])  # s-map coefficient
        coeff.m = colMeans(coeff[, -c(1, dim(coeff)[2]), drop = FALSE], na.rm = T)
        
        # S-map coefficients
        x1 = x
        x1[which(x1==0)] = NA
        x1[which(x1==1)] = coeff.m
        
        rho = block_lnlp_output$rho
        
        # test on the significance of rho
        n_pred = block_lnlp_output$num_pred
        t = rho*sqrt(n_pred-2)/sqrt(1-rho^2)  
        pvalue = 2*pt(-abs(t), df=n_pred-2)
        
        return(c(x1, theta.opt, rho, pvalue))  # c(coefficients, theta, rho, pvalue)
    })
    output = data.frame(t(coeff.mean))
    colnames(output) = c(names(data)[-1], 'theta', 'rho', 'pvalue')
    
    return(output)
}



