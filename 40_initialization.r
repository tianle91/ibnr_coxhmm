# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


empars_ini = function(shapes,
                      dat, 
                      dat_scalemodifier) {
    
    # shapes: c(shape1, shape2, shape3, ... shapem)
    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
                      
    nstates = length(shapes)
    nobs = length(dat)
    
    # initial probability is uniform
    mc_iniprob = rep(1/nstates,nstates)

    mc_transition = matrix(0, nrow = nstates, ncol = nstates)
    for (row in 1:nrow(mc_transition)) {
        prob = rep(1,nstates)
        # transition probability from self-state to self-state is equal to self-state to otherwise
        prob[row] = sum(prob) - prob[row]
        # regularize such that sumto1
        prob = prob / sum(prob)
        mc_transition[row,] = prob
    }
    
    # initial scale as proposed in Section 3.2.3
    scale = nstates / nobs
    scale = scale * sum(dat / dat_scalemodifier)
    scale = scale / sum(shapes)
    
    pars_list = list(shapes = shapes,
                     scale = scale)
    
    # return initial empars, which is a list of parameters needed for em to begin
    return(list(pars_list = pars_list,
                mc_iniprob = mc_iniprob,
                mc_transition = mc_transition))
}
# empars_ini(1:10, dat, dat_scalemodifier)


test = F
if (test) {
    set.seed(999471646)
    
    ndat = 1e2
    
    sample = c(0,1,2,3,4)
    dat = rep(sample, 1e2/length(sample))
    dat_scalemodifier = runif(ndat, min = .9, max = 1.1)
}