# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


single_shape_change = function(mc_em_output, 
                               comp_dex, change_by, 
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
    # RETURNS new em_output for a particular shape, changed by one step
    
    # denote best em_output
    best = mc_em_output

    # change original shape parameters
    empars_new = best$empars
    # cap new shape such that it is no smaller than abs(change_by)
    empars_new$pars_list$shapes[comp_dex] = max(abs(change_by), 
                                                empars_new$pars_list$shapes[comp_dex]+change_by)
                                                
    # get em_output for this new shape
    if (verbose) {
        cat('performing em on shape:', empars_new$pars_list$shapes, '\t')
    }
    new = mc_em(dat, dat_scalemodifier, empars_new, verbose = F)

    if (verbose) {
        cat('llike_old', best$loglike,
            'llike_new', new$loglike, '\t')
    }
    
    # return new iff new has higher loglike
    if (new$loglike > best$loglike) {
        if (verbose) { cat('new found\n') }
        return(new)
    } else {
        if (verbose) { cat('prev better\n') }
        return(best)
    }
}


single_shape_change_ttm = function(mc_em_output, 
                                   comp_dex, change_by, change_stepmax, 
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
    # RETURNS new em_output for a particular shape, changed until loglike dips
    
    best = mc_em_output
    temp = single_shape_change(best,
                               comp_dex, change_by, 
                               dat, 
                               dat_scalemodifier,
                               verbose = verbose)
                               
    # set up iterative
    step = 1
    while (!shapeisequal(temp, best) & 
           step <= change_stepmax) {
        # update best if step does not exceed max
        best = temp 
        # if step will not exceed max on next step, i.e. while loop still runs,
        # update new temp em_output
        if (step < change_stepmax) {
            step = step + 1
            temp = single_shape_change(best,
                                       comp_dex, change_by, 
                                       dat, 
                                       dat_scalemodifier,
                                       verbose = verbose)
        }
    }
    return(best)
}


shape_search_single = function(mc_em_output, 
                               stepsize = 1, 
                               stepmax = 10 * max(mc_em_output$empars$pars_list$shapes),
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
    # RETURNS new em_output after running increments and decrements as stated by leelin
    
    nstates = length(mc_em_output$empars$pars_list$shapes)

    # denote current best em_output
    best = mc_em_output
    changed = rep(F, nstates)
    
    for (k in nstates:1) {
        # test increments
        if (verbose) {cat('\ntest increment for comp #', k, '\n')}

        temp = single_shape_change_ttm(best, 
                                       comp_dex = k, 
                                       change_by = stepsize, change_stepmax = stepmax,
                                       dat, 
                                       dat_scalemodifier,
                                       verbose = verbose)
        
        # if temp has different shape than best, then new and better shape has been found
        if (!shapeisequal(temp, best)) {
            # replace best with this new em_output
            best = temp 
            changed[k] = T
        }
    }
    
    for (k in (1:nstates)[!changed]) {
        # only test decrement on non-incremented
        if (verbose) {cat('\ntest decrement for comp #', k, '\n')}
        
        temp = single_shape_change_ttm(best, 
                                       comp_dex = k, 
                                       change_by = -stepsize, change_stepmax = stepmax,
                                       dat,
                                       dat_scalemodifier,
                                       verbose = verbose)
                                       
        # if temp has different shape than best, then new and better shape has been found
        if (!shapeisequal(temp, best)) {
            # replace best with this new em_output
            best = temp 
        }
    }
    return(best)
}

shape_search = function(mc_em_output, 
                        stepsize = 1, 
                        stepmax = 10 * max(mc_em_output$empars$pars_list$shapes),
                        loopmax = 10^3,
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
    # RETURNS new em_output after running increments and decrements as stated by leelin
    
    # denote current best em_output
    best = mc_em_output
    
    for (loopdex in 1:loopmax) {
        if (verbose) {cat('\nshapesearch single for loopdex', loopdex, '\n')}
        temp = shape_search_single(best, 
                                   stepsize = stepsize, 
                                   stepmax = stepmax,
                                   dat = dat,
                                   dat_scalemodifier = dat_scalemodifier,
                                   verbose = verbose)
        if (!shapeisequal(temp, best)) {
            if (verbose) {cat('\nloopdex ', loopdex, 'shows improvement, run another loop\n')}
            # replace best with this new em_output
            best = temp 
        } else {
            if (verbose) {cat('\nloopdex ', loopdex, 'does not show improvement, return best\n')}
            return(best)
        }
    }
    if (verbose) {cat('\nloopmax of', loopmax, 'reached, return best\n')}
    return(best)
}