# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


loglike = function(dat, dat_scalemodifier, empars) {

    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    # empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    
    # make state probabilities, required to call logfwd function
    mc_genlogprob_output = mc_genlogprob_dnbinom(dat, dat_scalemodifier, empars$pars_list)
    # logfwd function also calculates log-likelihood
    mc_genlogfwd_output = mc_genlogfwd(mc_iniprob = empars$mc_iniprob, 
                                       mc_transition = empars$mc_transition, 
                                       mc_genlogprob_output)
    
    return(mc_genlogfwd_output$loglike)
}
# loglike(dat, dat_scalemodifier, empars_ini)

aic = function(dat, dat_scalemodifier, empars) {

    # dat: c(x1, x2, x3, ..., xn)
    # dat_scalemodifier: c(a1, a2, a3, ..., an)
    # empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    
    llike = loglike(dat, dat_scalemodifier, empars)
    
    nstates = length(empars$mc_iniprob)
    # number of freepars is (nstates^2 + nstates), agreed with dameng
    # but this should be computationally cheap and much more easily verfiable
    freepars = nstates +             # shapes
               1 +                   # scale 
               (nstates - 1) +       # ini_prob
               (nstates^2 - nstates) # transition matrix
    return(-2 * (llike - freepars))
}

test = F
if (test) {
    
    # 01 testing shit ------------------------------------------------------------------------------
    
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
    
    # 02 testing shit ------------------------------------------------------------------------------
    
    source('01_hmm_fwdbwd.r')
    
    empars_ini = list(pars_list = pars_list, mc_iniprob = mc_iniprob, mc_transition = mc_transition)
}