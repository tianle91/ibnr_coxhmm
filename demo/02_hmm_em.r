# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


mc_estep_statemarginal = function(mc_genlogfwd_output, 
                                  mc_genlogbwd_output) {
    
    # nrow = ndat, ncol = nstates
    ndat   = nrow(mc_genlogfwd_output$logfwd)
    nstate = ncol(mc_genlogfwd_output$logfwd)
    
    out = matrix(NA, nrow = ndat, ncol = nstate)
    for (j in 1:nstate) {
        # log version of the desired value is calculated
        # out[i, j] = z_{(t=i), j} in (3.11)
        value = mc_genlogfwd_output$logfwd[ , j] + 
                mc_genlogbwd_output$logbwd[j,  ] -
                mc_genlogfwd_output$loglike
        out[, j] = exp(value)
    }    
    return(out)
}
# mc_estep_statemarginal(mc_genlogfwd_output, mc_genlogbwd_output)


mc_mstep_next_hmmpars = function(mc_estep_statemarginal_output, 
                                 mc_genlogfwd_output, 
                                 mc_genlogbwd_output,
                                 mc_transition,
                                 mc_genlogprob_output) {
    
    nstate = ncol(mc_estep_statemarginal_output)
    ndat   = nrow(mc_estep_statemarginal_output)
    
    # follows (3.13)
    mc_iniprob_new = mc_estep_statemarginal_output[1,]
    
    # from (3.12, 3.14)
    mc_transition_new = matrix(NA, nrow = nstate, ncol = nstate)
    for (i in 1:nstate) {
        for (j in 1:nstate) {
            # this is log( (3.12) / gamma[i,j] )
            # here, we see log prob for neg-binom probs come in useful
            value = as.numeric(mc_genlogfwd_output$logfwd[1:(ndat-1), i     ]) + 
                    as.numeric(mc_genlogprob_output      [2:ndat    , j     ]) +
                    as.numeric(mc_genlogbwd_output$logbwd[j         , 2:ndat]) - 
                    mc_genlogfwd_output$loglike
            # take gamma[i,j] out and sum (3.12), then multiply in gamma[i,j]
            # obtain denominator of (3.14)
            value = sum(exp(value)) * mc_transition[i,j]
            mc_transition_new[i,j] = value
        }
    }
    # check for rowsums is 0
    rowsums_is0 = (1:nstate)[rowSums(mc_transition_new) == 0]
    if (length(rowsums_is0) > 0) {
        mc_transition_new[rowsums_is0, rowsums_is0] = diag(1, length(rowsums_is0))
    }
    # divby rowsums to get (3.14) proper
    mc_transition_new = mc_transition_new / rowSums(mc_transition_new)

    return(list(mc_iniprob = mc_iniprob_new,
                mc_transition = mc_transition_new))
}
# mc_estep_statemarginal_output = mc_estep_statemarginal(mc_genlogfwd_output, mc_genlogbwd_output)
# mc_estep_statejoint_output = mc_estep_statejoint(mc_genlogfwd_output, mc_genlogbwd_output, mc_transition, mc_genlogprob_output)
# mc_mstep_next_hmmpars(mc_estep_statemarginal_output, mc_estep_statejoint_output)


mc_mstep_next_scale = function(dat,
                               dat_scalemodifier,
                               pars_list,
                               mc_estep_statemarginal_output) {
    
    nstate = ncol(mc_estep_statemarginal_output)
    nobs   = nrow(mc_estep_statemarginal_output)

    # follows ztilde as stated in line above (3.16)
    ztilde = rowSums(mc_estep_statemarginal_output * matrix(rep(pars_list$shapes, nobs), 
                                                            byrow = T, nrow = nobs))

    # domain limits to be passed to uniroot
    # such that all terms in sum are positive
    scalemax = max(dat) / min(ztilde * dat_scalemodifier) 
    scalemax = min(1 / .Machine$double.eps, scalemax) # cap scalemax by some large value
    # cat('scalemax', scalemax, '\n')
    
    # uniroot objective function
    obj = function(scale){
        # equivalent expression from dameng, avoids the expression:
        # (ztilde a_t theta - n_t) in (3.16)
        numerator = sum(dat / (1 + dat_scalemodifier * scale))
        denominator = sum((ztilde * dat_scalemodifier) / (1 + dat_scalemodifier * scale))
        return(numerator / denominator - scale)
    }
    
    # finding root is very important here
    # interval reduces search range
    # request small tolerance 1e-32 
    # large maxiter of 1e6
    # stop and raise exception if not converged
    scale_new = uniroot(obj, interval = c(0, scalemax),
                        tol = .Machine$double.eps^2,
                        maxiter = 1e6,
                        check.conv = T)$root
    return(scale_new)
}
# mc_estep_statemarginal_output = mc_estep_statemarginal(mc_genlogfwd_output, mc_genlogbwd_output)
# mc_mstep_next_scale(dat, dat_scalemodifier, pars_list, mc_estep_statemarginal_output)


mc_emiter_single = function(dat,
                            dat_scalemodifier,
                            empars) {
                            
    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    # empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    
    pars_list     = empars$pars_list
    mc_iniprob    = empars$mc_iniprob
    mc_transition = empars$mc_transition
    
    # get probabilities necessary for E/M steps 
    mc_genlogprob_output = mc_genlogprob_dnbinom(dat, 
                                                 dat_scalemodifier, 
                                                 pars_list)
    mc_genlogfwd_output  = mc_genlogfwd(mc_iniprob, 
                                        mc_transition, 
                                        mc_genlogprob_output)
    mc_genlogbwd_output  = mc_genlogbwd(mc_iniprob, 
                                        mc_transition, 
                                        mc_genlogprob_output)
    
    # get estep results
    mc_estep_statemarginal_output = mc_estep_statemarginal(mc_genlogfwd_output, 
                                                           mc_genlogbwd_output)
    # get next parameters
    mc_hmm_pars_next              = mc_mstep_next_hmmpars(mc_estep_statemarginal_output, 
                                                          mc_genlogfwd_output, 
                                                          mc_genlogbwd_output,
                                                          mc_transition,
                                                          mc_genlogprob_output)
    mc_scale_next                 = mc_mstep_next_scale(dat, 
                                                        dat_scalemodifier,
                                                        pars_list, 
                                                        mc_estep_statemarginal_output)
    # check for nans in next iteration
    if (is.nan(sum(c(mc_scale_next,
                     as.numeric(mc_hmm_pars_next$mc_iniprob),
                     as.numeric(mc_hmm_pars_next$mc_transition))))) {
        print('called empars')
        print(empars)
        print('new scale')
        print(mc_scale_next)
        print('new hmm pars')
        print(mc_hmm_pars_next)
        stop('nans in next iteration\n', call. = T)
    }
    
    # denote new shits
    pars_list_new       = pars_list
    pars_list_new$scale = mc_scale_next
    mc_iniprob_new      = mc_hmm_pars_next$mc_iniprob
    mc_transition_new   = mc_hmm_pars_next$mc_transition
    
    # return new empars
    return(list(pars_list = pars_list_new,
                mc_iniprob = mc_iniprob_new,
                mc_transition = mc_transition_new))
}
# empars_ini = list(pars_list = pars_list, mc_iniprob = mc_iniprob, mc_transition = mc_transition)
# empars_new = mc_emiter_single(dat, dat_scalemodifier, empars_ini)


mc_diff = function(empars_prev, empars_new) {

    # empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    
    # scale 1norm
    scale_diff = empars_prev$pars_list$scale - empars_new$pars_list$scale
    
    # iniprob 1norm
    mc_iniprob_diff = empars_prev$mc_iniprob - empars_new$mc_iniprob
    mc_iniprob_diff = as.numeric(mc_iniprob_diff)
    
    # transition matrix 1norm
    mc_transition_diff = empars_prev$mc_transition - empars_new$mc_transition
    mc_transition_diff = as.numeric(mc_transition_diff)
    
    return(sum(abs(c(scale_diff, 
                     mc_iniprob_diff,
                     mc_transition_diff))))
}
# empars_ini = list(pars_list = pars_list, mc_iniprob = mc_iniprob, mc_transition = mc_transition)
# empars_new = mc_emiter_single(dat, dat_scalemodifier, empars_ini)
# mc_diff(empars_new$empars, empars_ini)


mc_em = function(dat,
                 dat_scalemodifier, 
                 empars,
                 ntrials = 1e6,
                 reltol = 1e-6, 
                 verbose = T) {
                 
    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    # empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    
    # set up current best
    empars_current = empars
    loglike_current = loglike(dat, dat_scalemodifier, empars_current)
    
    # set up stopping condition counters
    trialdex = 1
    reldiff = Inf
    
    while (trialdex <= ntrials) {
        
        # get outputs
        empars_new = mc_emiter_single(dat, dat_scalemodifier, empars_current)
        loglike_new = loglike(dat, dat_scalemodifier, empars_new)
        
        # absolute difference 
        reldiff = mc_diff(empars_new, empars_current)
        
        if (verbose) {
            cat('iter', trialdex,
                'llk', loglike_current,
                'rdif', reldiff,
                '\n')
        }
        
        # stop if loglike dips
        # give allowance of 6th decimal point difference
        # saw case of 12th dp difference triggering, shouldn't be a huge problem
        if (loglike_new < (loglike_current - 1e-6)) {
            options(digits = 22)
            print('loglike_current')
            print(loglike_current)
            print('loglike_new')
            print(loglike_new)
            print('empars_current')
            print(empars_current)
            stop('loglike dipped', call. = TRUE)
        }

        # stop if reltol reached
        if (reldiff < reltol) {
            if (verbose) {
                cat('reldiff of', reltol, 'reached in', trialdex, 'iterations.\n')
            }
            return(list(empars = empars_new,
                        loglike = loglike_new,
                        aic = aic(dat, dat_scalemodifier, empars_new)))
        }

        # set up next iteration
        trialdex = trialdex + 1
        empars_current = empars_new
        loglike_current = loglike_new
    }
    
    # stop if reached max trials
    if (verbose) {
        cat('ntrials of', ntrials,
            'reached, with finalreldiff of', reldiff,
            'loglike', loglike_current,
            '\n')
    }
    return(list(empars = empars_current,
                loglike = loglike_current,
                aic = aic(dat, dat_scalemodifier, empars_current)))
}
# empars_ini = list(pars_list = pars_list, mc_iniprob = mc_iniprob, mc_transition = mc_transition)
# mc_em(dat, dat_scalemodifier, empars_ini, verbose = T)


test = F
if (test) {
    
    # 01 testing shit ------------------------------------------------------------------------------
    
    set.seed(999471646)
    
    dat = rnbinom(10, size = 5, prob = .5)
    dat_scalemodifier = 1 + rnbinom(10, size = 1, prob = .3)
    
    mc_iniprob = runif(5)
    mc_iniprob = mc_iniprob / sum(mc_iniprob)
    # mc_iniprob = c(1,0)
    
    mc_transition = matrix(0, nrow = 5, ncol = 5)
    for (row in 1:nrow(mc_transition)) {
        prob = runif(5)
        prob = prob / sum(prob)
        mc_transition[row,] = prob
    }
    # mc_transition = rbind(c(0,1),
    #                       c(0,1))
    rm(row,prob)
    
    shapes = 1:5
    # shapes = c(1,1)
    scale = 1
    pars_list = list(shapes = shapes,
                     scale = scale)
    rm(shapes, scale)
    
    
    # 02 testing shit ------------------------------------------------------------------------------
    
    source('01_hmm_fwdbwd.r')
    source('01a_hmm_loglike.r')
    
    mc_genlogprob_output = mc_genlogprob_dnbinom(dat, dat_scalemodifier, pars_list)
    mc_genlogfwd_output  = mc_genlogfwd(mc_iniprob, mc_transition, mc_genlogprob_output)
    mc_genlogbwd_output  = mc_genlogbwd(mc_iniprob, mc_transition, mc_genlogprob_output)
}