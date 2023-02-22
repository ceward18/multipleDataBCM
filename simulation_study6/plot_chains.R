
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

plot_alarm_chains <- function(chains, alarmString) {
    yAlarm1 <- chains[[1]][,grep(alarmString, colnames(chains[[1]]))]
    yAlarm2 <- chains[[2]][,grep(alarmString, colnames(chains[[2]]))]
    yAlarm3 <- chains[[3]][,grep(alarmString, colnames(chains[[3]]))]
    
    yAlarmPost <- cbind.data.frame(x = rep(1:ncol(yAlarm1), 3),
                                   Chain = rep(1:3, each = ncol(yAlarm1)),
                                   mean = c(colMeans(yAlarm1),
                                            colMeans(yAlarm2),
                                            colMeans(yAlarm3)),
                                   lower = c(apply(yAlarm1, 2, quantile, probs = 0.025),
                                             apply(yAlarm2, 2, quantile, probs = 0.025),
                                             apply(yAlarm3, 2, quantile, probs = 0.025)),
                                   upper = c(apply(yAlarm1, 2, quantile, probs = 0.975),
                                             apply(yAlarm2, 2, quantile, probs = 0.975),
                                             apply(yAlarm3, 2, quantile, probs = 0.975)))
    
    yAlarmPost$Chain <- factor(yAlarmPost$Chain)
    
    ggplot(yAlarmPost, aes(x = x, y = mean,
                           ymin = lower, ymax = upper,
                           group = Chain, col = Chain, fill = Chain)) +
        geom_line(size = 1) + 
        geom_ribbon(alpha = 0.2) +
        theme_bw() +
        ggtitle('Alarm function') + 
        scale_color_manual(values = c('black', 'red', 'blue')) + 
        scale_fill_manual(values = c('black', 'red', 'blue')) +
        labs(x = '',  y = 'Alarm') +
        theme(plot.title = element_text(h = 0.5, size = 14),
              legend.title = element_text(h = 0.5, size = 13),
              legend.text = element_text(h = 0.5, size = 12),
              axis.text.y = element_text(h = 0.5, size = 12),
              axis.title.y = element_text(h = 0.5, size = 13),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())
    
}

plot_r0_chains <- function(chains) {
    yAlarm1 <- chains[[1]][,grep('R0', colnames(chains[[1]]))]
    yAlarm2 <- chains[[2]][,grep('R0', colnames(chains[[2]]))]
    yAlarm3 <- chains[[3]][,grep('R0', colnames(chains[[3]]))]
    
    yAlarmPost <- cbind.data.frame(x = rep(1:ncol(yAlarm1), 3),
                                   Chain = rep(1:3, each = ncol(yAlarm1)),
                                   mean = c(colMeans(yAlarm1),
                                            colMeans(yAlarm2),
                                            colMeans(yAlarm3)),
                                   lower = c(apply(yAlarm1, 2, quantile, probs = 0.025),
                                             apply(yAlarm2, 2, quantile, probs = 0.025),
                                             apply(yAlarm3, 2, quantile, probs = 0.025)),
                                   upper = c(apply(yAlarm1, 2, quantile, probs = 0.975),
                                             apply(yAlarm2, 2, quantile, probs = 0.975),
                                             apply(yAlarm3, 2, quantile, probs = 0.975)))
    
    yAlarmPost$Chain <- factor(yAlarmPost$Chain)
    
    ggplot(yAlarmPost, aes(x = x, y = mean,
                           ymin = lower, ymax = upper,
                           group = Chain, col = Chain, fill = Chain)) +
        geom_line(size = 1) + 
        geom_ribbon(alpha = 0.2) +
        theme_bw() +
        ggtitle('R0') + 
        scale_color_manual(values = c('black', 'red', 'blue')) + 
        scale_fill_manual(values = c('black', 'red', 'blue')) +
        labs(x = 'Epidemic Time',  y = 'R0(t)') +
        theme(plot.title = element_text(h = 0.5, size = 14),
              legend.title = element_text(h = 0.5, size = 13),
              legend.text = element_text(h = 0.5, size = 12),
              axis.text.y = element_text(h = 0.5, size = 12),
              axis.title.y = element_text(h = 0.5, size = 13),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())
    
}

plot_simple_chains <- function(chains, modelType) {
    
    par(mfrow = c(2,5))
    
    plot_param_chains(chains, 'probDetect')
    abline(h = 0.25, lty = 2, lwd = 2)
    plot_param_chains(chains, 'beta')
    abline(h = 0.7, lty = 2, lwd = 2)
    plot_param_chains(chains, 'delta[1]')
    if (modelType == 'inc') {
        abline(h = 0.6, lty = 2, lwd = 2)
    } else if (modelType == 'death') {
        abline(h = 0.1, lty = 2, lwd = 2)
    } else if (modelType == 'equal') {
        abline(h = 0.35, lty = 2, lwd = 2)
    }
    plot_param_chains(chains, 'delta[2]')
    if (modelType == 'inc') {
        abline(h = 0.1, lty = 2, lwd = 2)
    } else if (modelType == 'death') {
        abline(h = 0.6, lty = 2, lwd = 2)
    } else if (modelType == 'equal') {
        abline(h = 0.35, lty = 2, lwd = 2)
    }
    plot_param_chains(chains, 'x0C')
    abline(h = 100, lty = 2, lwd = 2)
    plot_param_chains(chains, 'x0D')
    abline(h = 15, lty = 2, lwd = 2)
    plot_param_chains(chains, 'nuC')
    abline(h = 3, lty = 2, lwd = 2)
    plot_param_chains(chains, 'nuD')
    abline(h = 3, lty = 2, lwd = 2)
    plot_param_chains(chains, 'k')
    plot_param_chains(chains, 'w0')
    
    
}

plot_inc_chains <- function(chains, modelType) {
    
    par(mfrow = c(2,4))
    
    plot_param_chains(chains, 'probDetect')
    abline(h = 0.25, lty = 2, lwd = 2)
    plot_param_chains(chains, 'beta')
    abline(h = 0.7, lty = 2, lwd = 2)
    plot_param_chains(chains, 'deltaC')
    if (modelType == 'inc') {
        abline(h = 0.6, lty = 2, lwd = 2)
    } else if (modelType == 'death') {
        abline(h = 0.1, lty = 2, lwd = 2)
    } else if (modelType == 'equal') {
        abline(h = 0.35, lty = 2, lwd = 2)
    }
    plot_param_chains(chains, 'x0C')
    abline(h = 100, lty = 2, lwd = 2)
    plot_param_chains(chains, 'nuC')
    abline(h = 3, lty = 2, lwd = 2)
    plot_param_chains(chains, 'k')
    plot_param_chains(chains, 'w0')
    
    
}

plot_full_chains <- function(chains, modelType) {
    
    par(mfrow = c(3,4))
    
    plot_param_chains(chains, 'probDetect')
    abline(h = 0.25, lty = 2, lwd = 2)
    plot_param_chains(chains, 'beta')
    abline(h = 0.7, lty = 2, lwd = 2)
    plot_param_chains(chains, 'delta[1]')
    if (modelType == 'inc') {
        abline(h = 0.6, lty = 2, lwd = 2)
    } else if (modelType == 'death') {
        abline(h = 0.1, lty = 2, lwd = 2)
    } else if (modelType == 'equal') {
        abline(h = 0.35, lty = 2, lwd = 2)
    }
    plot_param_chains(chains, 'delta[2]')
    if (modelType == 'inc') {
        abline(h = 0.1, lty = 2, lwd = 2)
    } else if (modelType == 'death') {
        abline(h = 0.6, lty = 2, lwd = 2)
    } else if (modelType == 'equal') {
        abline(h = 0.35, lty = 2, lwd = 2)
    }
    plot_param_chains(chains, 'x0C')
    abline(h = 100, lty = 2, lwd = 2)
    plot_param_chains(chains, 'x0D')
    abline(h = 15, lty = 2, lwd = 2)
    plot_param_chains(chains, 'nuC')
    abline(h = 3, lty = 2, lwd = 2)
    plot_param_chains(chains, 'nuD')
    abline(h = 3, lty = 2, lwd = 2)
    plot_param_chains(chains, 'gamma1')
    abline(h = 0.2, lty = 2, lwd = 2)
    plot_param_chains(chains, 'gamma2')
    abline(h = 0.2, lty = 2, lwd = 2)
    plot_param_chains(chains, 'lambda')
    abline(h = 0.1, lty = 2, lwd = 2)
    plot_param_chains(chains, 'phi')
    abline(h = 0.1, lty = 2, lwd = 2)
    
    
}



plot_simple_chains(resThree, dataType_i)

plot_inc_chains(resThree, dataType_i)

plot_full_chains(resThree, dataType_i)