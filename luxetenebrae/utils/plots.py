import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

def set_plot_style():
    """
    Set the global plot style for matplotlib.
    This function configures the default aesthetics for plots.
    """
    plt.rcParams['axes.grid'] = False  # Disable grid lines
    plt.rcParams['grid.alpha'] = 0.5  # Set grid transparency
    plt.rcParams['savefig.dpi'] = 300  # High resolution for saved figures
        # Set global font to Times New Roman, 11pt
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman"]
    plt.rcParams["font.size"] = 11
    # Set default colormap to plasma
    plt.rcParams["image.cmap"] = "plasma"


# Formatter for 2 significant figures
def two_sig_figs(val, pos=None):
    if val == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(val))))
    coeff = round(val / 10**exponent, 1)
    return r"${} \times 10^{{{}}}$".format(coeff, exponent)

# Example: plot stellar type evolution (from Stage_03a_generate_plots.py)
def plot_time_evolution(array_x, array_y, y_labels,legend_labels, path_to_save, fill=None, x_labels=None, idx=None, info=None):
    rows = len(array_y)
    cols = len(array_x)
    fig, ax = plt.subplots(rows, cols, sharex=True, figsize=(3.5, cols*3.5))
    for i in range(rows*cols):
        for j in range(len(array_y[i])):
            x = array_x[j]
            y = array_y[i][j]
            label = legend_labels[j]
            ax[i].plot(x, y, label=label)
        ax[i].xaxis.set_major_formatter(FuncFormatter(two_sig_figs))
        ax[i].yaxis.set_major_formatter(FuncFormatter(two_sig_figs))
        ax[i].legend()
        ax[i].set_ylabel(y_labels[i])
        if x_labels is not None:
            ax[i].set_xlabel(x_labels[i])



    if x_labels is None:
        # If no x_labels provided, default to 'Time [Myr]'
        plt.set_xlabel('Time [Myr]')

    if info:
        for i, text in enumerate(info):
            ax[i].text(0.05, 0.95, text, transform=ax[i].transAxes, verticalalignment='top')

    if idx is not None:
        filename = f'{path_to_save}/stellar_type_evolution_{idx}.png'
    else:
        filename = f'{path_to_save}/stellar_type_evolution.png'
    plt.savefig(filename, bbox_inches='tight')
    plt.close()

# Example: plot histogram of properties (from Stage_03b_detailed_evolution.py)
def plot_bhms_properties_histogram(rows, cols, array_x, array_y, y_labels,legend_labels, path_to_save, x_labels=None, x_scale=None, y_scale=None, idx=None, info=None):
    fig, ax = plt.subplots(rows, cols, figsize=(3.5*rows, 3.5*cols))
    for i in range(rows*cols):
        x = array_x[i]
        y = array_y[i]
        label = legend_labels[i]
        ax[i].hist(y, bins=20, color='coral')
        ax[i].set_title(label)
        ax[i].xaxis.set_major_formatter(FuncFormatter(two_sig_figs))
        ax[i].yaxis.set_major_formatter(FuncFormatter(two_sig_figs))
    plt.tight_layout()
    plt.savefig(path_to_save + 'bhms_properties_histogram.png')
    plt.close(fig)

# Add more plotting functions as needed, following this pattern.
