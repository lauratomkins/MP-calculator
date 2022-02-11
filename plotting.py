import numpy as np
import matplotlib.pyplot as plt


def plot_curves(ax, model_times, plot_dict, use_first=True, first_val=None, ylab=None, xlab=None,
                title=None, yscale=None, xscale=None):

    for key in plot_dict:
        plot_var = plot_dict[key]['data']
        if len(plot_var) < len(model_times):
            if use_first:
                plot_var = np.insert(plot_var, 0, plot_var[0])
            else:
                plot_var = np.insert(plot_var, 0, first_val)
        plot_var = np.ma.masked_values(plot_var, 0)
        ax.plot(model_times, plot_var, label=plot_dict[key]['label'], linewidth=2)

    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    #ax1.set_ylim([0,60])
    if yscale == 'log': ax.set_yscale('log')
    if xscale == 'log': ax.set_xscale('log')
    ax.grid()
    ax.set_title(title, loc='right')
    ax.legend()