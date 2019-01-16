# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


em_compno_search = function(mc_em_output, 
                            dat,
                            dat_scalemodifier, 
                            shapestep = 1,
                            min_comps = 1, 
                            max_loop = length(mc_em_output$empars$pars_list$shapes),
                            verbose = T) {
                            
    # mc_em_output = list(empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition),
    #                     loglike, aic)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    #
    # RETURNS find em_output with best and min number of components by aic
    
    # denote current best
    best = mc_em_output

    # do shapesearch for current best
    if (verbose) {
        cat('getting best shape for current shape:', 
            best$empars$pars_list$shapes, '\n')
    }
    best = shape_search(best, 
                        stepsize = shapestep,
                        dat = dat, dat_scalemodifier = dat_scalemodifier, 
                        verbose = verbose)
    nstates_best = length(best$empars$pars_list$shapes)
    
    # loop through and evaluate competitors
    for (loopdex in 1:max_loop) {
        
        # denote competitor by attempting to reduce comps
        if (verbose) {
            cat('\ndrop min longrun proportion component from best and find best\n')
        }
        temp = reduce_comps(best, 
                            dat, dat_scalemodifier,
                            shapestep = shapestep,
                            verbose = verbose)
        nstates_temp = length(temp$empars$pars_list$shapes)
        
        if (verbose) {
            cat('loopdex:', loopdex,
                'best nstates:', nstates_best, 
                'temp nstates:', nstates_temp, 
                '\n')
        }
        
        if (nstates_temp < nstates_best) {
            best = temp
            nstates_best = nstates_temp
        } else {
            if (verbose) {
                cat('min_comps reached, best nstates:', nstates_best,
                    '\n')
            }
            return(best)
        }
    }
    
    # stop if max loop reached
    if (verbose) {
        cat('maxloop finished, best nstates:', nstates_best,
            '\n')
    }
    return (best)   
}