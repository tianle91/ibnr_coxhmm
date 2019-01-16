hmm_sample = function(maxt, iniprob, transprob, shapes, scale, scalemod=NULL) {

    # generates samples and latent state for a nbinom-hmm
    #
    # Args:
    #     maxt: number of time steps
    #     iniprob: vector of initial probabilities
    #     transprob: matrix of markov transition probabilities
    #     shapes: vector of shapes corresponding to emission 
    #     scale: common scale for emission parameters
    #     scalemod: vector of scale modifiers. if NULL, this is vector of 1s.
    #
    # Returns:
    #     list of samples and latent states

    if (is.null(scalemod)) {
        scalemod = rep(1, maxt)
    } else if (not(length(scalemod) == maxt)) {
    	stop('length of scalemod provided is not maxt!')
    }
    
    nstates = length(shapes)
    x = rep(0, maxt)
    c = rep(0, maxt)
    
    tdex = 1

    while (tdex <= maxt) {

    	# sample latent state
        if (tdex == 1) {
            statenow = sample(1:nstates, 1, prob = iniprob)
        } else {
            statenow = sample(1:nstates, 1, prob = transprob[statenow,])
        }
        c[tdex] = statenow

        # get current emission parameters
        shapenow = shapes[statenow]
        scalenow = scalemod[tdex]*scale
        probnow = 1 / (1 + scalenow)

        # random sample from emission parameters
        x[tdex] = rnbinom(1, size = shapenow, prob = probnow)
        tdex = tdex + 1
    }

    return (list('dat'=x, 'latent'=c))
}