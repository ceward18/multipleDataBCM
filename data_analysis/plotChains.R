
plot_param_chains <- function(chains, param) {
    ylims <- c(min(chains[[1]][,param],
                   chains[[2]][,param],
                   chains[[3]][,param]),
               max(chains[[1]][,param],
                   chains[[2]][,param],
                   chains[[3]][,param]))
    
    plot(chains[[1]][,param], type = 'l', ylim = ylims, main = param, col = 'grey35')
    lines(chains[[2]][,param], col = 'palevioletred')
    lines(chains[[3]][,param], col = 'deepskyblue')
}

plot_simple_chains <- function(chains) {
    
    par(mfrow = c(3,4))
    
    plot_param_chains(chains, 'beta')
    plot_param_chains(chains, 'delta[1]')
    plot_param_chains(chains, 'delta[2]')
    plot_param_chains(chains, 'x0C')
    plot_param_chains(chains, 'x0D')
    plot_param_chains(chains, 'nuC')
    plot_param_chains(chains, 'nuD')
    plot_param_chains(chains, 'k')
    plot_param_chains(chains, 'w0')
    plot_param_chains(chains, 'comp_init[1]')
    plot_param_chains(chains, 'comp_init[2]')
    plot_param_chains(chains, 'comp_init[3]')
    
    
}

plot_full_chains <- function(chains) {
    
    par(mfrow = c(4,4))
    
    plot_param_chains(chains, 'beta')
    plot_param_chains(chains, 'delta[1]')
    plot_param_chains(chains, 'delta[2]')
    plot_param_chains(chains, 'x0C')
    plot_param_chains(chains, 'x0D')
    plot_param_chains(chains, 'nuC')
    plot_param_chains(chains, 'nuD')
    plot_param_chains(chains, 'comp_init[1]')
    plot_param_chains(chains, 'comp_init[2]')
    plot_param_chains(chains, 'comp_init[3]')
    plot_param_chains(chains, 'comp_init[4]')
    plot_param_chains(chains, 'comp_init[5]')
    plot_param_chains(chains, 'gamma1')
    plot_param_chains(chains, 'gamma2')
    plot_param_chains(chains, 'lambda')
    plot_param_chains(chains, 'phi')
    
    
}

plot_inc_chains <- function(chains) {
    
    par(mfrow = c(3,3))
    
    plot_param_chains(chains, 'beta')
    plot_param_chains(chains, 'deltaC')
    plot_param_chains(chains, 'x0C')
    plot_param_chains(chains, 'nuC')
    plot_param_chains(chains, 'k')
    plot_param_chains(chains, 'w0')
    plot_param_chains(chains, 'comp_init[1]')
    plot_param_chains(chains, 'comp_init[2]')
    plot_param_chains(chains, 'comp_init[3]')
    
    
}


plot_full_chains(resThree)





