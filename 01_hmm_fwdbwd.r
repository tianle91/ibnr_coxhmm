# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


mc_genlogprob_dnbinom = function(dat, 
                                 dat_scalemodifier,
                                 pars_list,
                                 verbose = F) {
                                 
    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    # pars_list = list(shapes, scale)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale

    nvals   = length(dat)
    nstates = length(pars_list$shapes)
    
    shapes = pars_list$shapes
    scales = pars_list$scale * dat_scalemodifier

    # generate matrix with columns corresponding to pars_list    
    out = matrix(0, nrow = nvals, ncol = nstates)
    for (datdex in 1:nvals) {
        # note that out[i, j] = logprob(xi | state = j)
        
        # alternative parameterization follows that derived from discretized gamma
        prob = 1 / (1 + scales[datdex])
        # out_thisrow = logprob(xi | state = {1, 2, 3, ... })
        # note that log allows for later E/M step calculations to not underflow often
        out_thisrow = dnbinom(dat[datdex], size = shapes, p = prob, log = T)
        
        # report modes given parameter, easy for fitting trivial hmms for max-loglike, 
        # where obs jump from one mode to another predictably
        if (verbose) {
            for (size in shapes) {
                cat('prob', prob,
                    'size', size,
                    'mode', floor(prob * (size - 1) / (1 - prob)),
                    '\n')
            }
        }
        out[datdex, ] = out_thisrow
    }
    return(out)
}
# mc_genlogprob_dnbinom(dat, dat_scalemodifier, pars_list)


mc_genlogfwd = function(mc_iniprob,
                        mc_transition,
                        mc_genlogprob_output) {
    
    nvals   = nrow(mc_genlogprob_output)
    nstates = ncol(mc_genlogprob_output)
    
    # preassign desired log(alpha) matrix
    logfwd_m = matrix(0, nrow = nvals, ncol = nstates)
    # preassign scaling vector
    logfwdsums = numeric(nvals)
    
    # initial assignments
    # note that probs given state is in log
    fwd       = t(mc_iniprob) * exp(mc_genlogprob_output[1,]) 
    logfwdsum = log(sum(fwd))
    scaled    = fwd / exp(logfwdsum)
    
    # write initial assignments
    logfwd_m[1,] = log(fwd)
    logfwdsums[1] = logfwdsum
    
    # recursive portion
    for (row in 2:nvals) {
    
        # useful row vector
        m = scaled %*% mc_transition * exp(mc_genlogprob_output[row,])
        
        # get new versions
        logfwdsum_new = logfwdsum + log(sum(m))
        scaled_new = exp(logfwdsum - logfwdsum_new) * m
        logfwd_new = log(scaled_new) + rep(logfwdsum_new, nstates)
        
        # write new
        logfwd_m[row,] = logfwd_new
        logfwdsums[row] = logfwdsum_new
        
        # update initial assignments
        scaled = scaled_new
        logfwdsum = logfwdsum_new
    }
    
    # get loglike (from dameng) with scaling modification
    scalingconst = max(logfwd_m[nvals,])
    loglike      = scalingconst + log(sum(exp(logfwd_m[nvals,] - scalingconst)))
    
    return (list(logfwd = logfwd_m,
                 logfwdsums = logfwdsums,
                 loglike = loglike))
}
# mc_genlogprob_output = mc_genlogprob_dnbinom(dat, dat_scalemodifier, pars_list)
# mc_genlogfwd(mc_iniprob, mc_transition, mc_genlogprob_output)


mc_genlogbwd = function(mc_iniprob,
                        mc_transition,
                        mc_genlogprob_output) {
    
    nvals = nrow(mc_genlogprob_output)
    nstates = ncol(mc_genlogprob_output)
    
    # preassign desired log(beta) matrix
    logbwd_m = matrix(0, nrow = nstates, ncol = nvals)
    # preassign scaling vector
    logbwdsums = numeric(nvals)
    
    # initial assignments
    # note that we start from the last column
    bwd       = matrix(rep(1,nstates), ncol = 1)
    logbwdsum = log(sum(bwd))
    scaled    = bwd / exp(logbwdsum)
    
    # write initial assignments
    logbwd_m[, nvals] = log(bwd)
    logbwdsums[nvals] = logbwdsum
    
    # recursive portion
    for (col in (nvals-1):1) {
    
        # useful col vector
        m = mc_transition %*% diag(exp(mc_genlogprob_output[col + 1,])) %*% scaled
        
        # get new versions
        logbwdsum_new = logbwdsum + log(sum(m))
        scaled_new = exp(logbwdsum - logbwdsum_new) * m
        logbwd_new = log(scaled_new) + matrix(rep(logbwdsum_new, nstates), ncol = 1)
        
        # write new
        logbwd_m[, col] = logbwd_new
        logbwdsums[col] = logbwdsum_new
        
        # update initial assignments
        scaled = scaled_new
        logbwdsum = logbwdsum_new
    }
    return (list(logbwd = logbwd_m,
                 logbwdsums = logbwdsums))
}
# mc_genlogprob_output = mc_genlogprob_dnbinom(dat, dat_scalemodifier, pars_list)
# mc_genlogbwd(mc_iniprob, mc_transition, mc_genlogprob_output)


test = F
if (test) {
    set.seed(999471646)
    
    dat = rnbinom(10, size = 5, prob = .5)
    dat_scalemodifier = 1 + rnbinom(10, size = 1, prob = .3)
    
    mc_iniprob = runif(5)
    mc_iniprob = mc_iniprob / sum(mc_iniprob)
    
    mc_transition = matrix(0, nrow = 5, ncol = 5)
    for (row in 1:nrow(mc_transition)) {
        prob = runif(5)
        prob = prob / sum(prob)
        mc_transition[row,] = prob
    }
    rm(row,prob)
    
    shapes = 1:5
    scale = 1
    pars_list = list(shapes = shapes,
                     scale = scale)
    rm(shapes, scale)
}