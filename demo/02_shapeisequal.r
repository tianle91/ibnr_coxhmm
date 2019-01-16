# --------------------------------------------------------------------------------------------------
# CHEN, TIANLE
# tianle.chen@mail.utoronto.ca
# --------------------------------------------------------------------------------------------------


shapeisequal = function(mc_em_output1, mc_em_output2) {

    # mc_em_output = list(empars = list(pars_list = list(shapes, scale), mc_iniprob, mc_transition),
    #                     loglike, aic)
    #   shapes: c(shape1, shape2, shape3, ... shapem)
    #   scale: scale
    
    shape1 = mc_em_output1$empars$pars_list$shapes
    shape2 = mc_em_output2$empars$pars_list$shapes
    
    # return TRUE iff all items are equal
    if (length(shape1) != length(shape2)) {
        return (F)
    } else {
        return (prod(shape1==shape2)==1)
    }
}