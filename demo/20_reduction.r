# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


drop_minweight = function(mc_em_output, 
                          dat, 
                          dat_scalemodifier,
                          verbose = F) {
    
    # mc_em_output = list(empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition),
    #                     loglike, aic)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    #
    # RETURNS new em_output with min longrun probability state dropped, if has stationary states
    
    # denote current best em_output
    best = mc_em_output
    
    # denote useful variables
    nstates = length(best$empars$pars_list$shapes)    
    mc_transition = best$empars$mc_transition
    
	# from dameng
	inv_m = t(diag(nstates) - mc_transition + 1)
	if (abs(det(inv_m)) <= .Machine$double.eps) {
	    # not invertible, thus use less desirable alternative
		stationary = matrix(best$empars$mc_iniprob, nrow = 1)
		for (i in 1:1e6) {
		    stationary = stationary %*% best$empars$mc_transition
		}
	} else {
	    # actual stationary distribution
	    stationary = solve(inv_m, rep(1, nstates))
	}
	
    # identify minimum proportion for longrun probability stat
    minweight = min(stationary)
    kickdex = (1:nstates)[stationary == minweight]
    keep    = (1:nstates)[stationary != minweight]
    if (length(kickdex) > 1) {
        keep = sort(c(kickdex[2:length(kickdex)], keep))
        kickdex = kickdex[1]
    }
    
    if (verbose) {
        cat('minweight:', minweight,
            'kickdex:', kickdex, 
            'nstates:', nstates,
            '\n')
        print(keep)
        print(stationary)
    }
    
    if (length(keep) > 1) {
        # if number of states to keep is more than 1
        # get newshape 
        newshape = best$empars$pars_list$shapes[keep]
        # note that initial probbilities and transitional probabilities are re-initialized
        # and do not carry over from previous estimates (because messy)
        empars   = empars_ini(newshape, dat, dat_scalemodifier)
        out      = mc_em(dat, dat_scalemodifier, empars, verbose = F)
        return (out)
    } else {
        # otherwise, do not reduce shapes and return current best
        return (best)
    }
    
}
# drop_minweight(mc_em_output, dat, dat_scalemodifier, verbose = T)


reduce_comps = function(mc_em_output, 
                        dat,
                        dat_scalemodifier, 
                        shapestep = 1,
                        verbose = F) {
                        
    # mc_em_output = list(empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition),
    #                     loglike, aic)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    #
    # RETURNS find best em_output with reduced number of states with shapesearch done
    # and compare aics to give best em_output for previous number of states and new reduced number
    
    best = mc_em_output

    # get new em_output from reducing number of states
    new = drop_minweight(best,
                         dat, 
                         dat_scalemodifier,
                         verbose = verbose)
    
    if (shapeisequal(new, best)) {
        return (best)
    }    
    # do shape_search for this new one
    new = shape_search(new,
                       stepsize = shapestep,
                       dat = dat, 
                       dat_scalemodifier = dat_scalemodifier,
                       verbose = verbose)

    # evaluate this new shit    
    improved = new$aic < best$aic
    if (verbose) {
        cat('old aic:', best$aic,
            'new aic:', new$aic,
            'improved:', improved,
            '\n')
    }
    
    if (improved) {
        return(new)
    } else {
        return(best)
    }
}