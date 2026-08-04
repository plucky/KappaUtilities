# Walter Fontana, 2026

import csv
import seaborn as sns
import matplotlib.pyplot as plt

# plt.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['font.family'] = 'Helvetica'
# mpl.rcParams['font.size'] = 18
# Avoid mathtext (renders as embedded Type-3 paths, not real text)
# for tick labels, and avoid the unicode minus sign, which isn't in
# the core-font encoding.
plt.rcParams['axes.formatter.use_mathtext'] = False
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['font.size'] = 18
# plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams["figure.constrained_layout.h_pad"] = 0
plt.rcParams["figure.constrained_layout.w_pad"] = 0
plt.rcParams["figure.constrained_layout.hspace"] = 1
plt.rcParams["figure.constrained_layout.wspace"] = 0
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True


def load_csv(file):
    x = []
    y = []
    with open(file) as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        x_name, y_name = reader.fieldnames
        for row in reader:
            x.append(float(row[x_name]))
            y.append(float(row[y_name]))
    return x, y


def line_plot(data=None, X=None, Y=None, hue=None, loop=None, size=(10, 5), ylim=None, xlim=None, legend=False,
              xlog=False, ylog=False, color='b', marker='o', msize=6, mcolor='red', drawstyle=None, linewidth=2,
              ax=None, plotname=None):
    #
    # don't set this if using PyCharm plot settings
    # sns.set_theme()
    # sns.set_style("whitegrid", {'axes.edgecolor': 'black', 'axes.linewidth': 1, 'grid.linestyle': ':'})
    # sns.set_style({'axes.linewidth': 1, 'grid.linestyle': ':'})

    if not ax:
        for i in plt.get_fignums():
            plt.close(plt.figure(i))
        fig, ax = plt.subplots(1, 1, figsize=size)
        # fig.patch.set_alpha(0.)
    else:
        # ax.patch.set_alpha(0.)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.)
        # ax.set_xlim(xlim)
        # ax.set_ylim(ylim)
        # ax.grid(True, alpha=0)

    colors = ['r', 'g', 'b', 'c', 'm', 'k', 'brown', 'darkorange', 'darkturquoise', 'slategray']

    if not loop:
        sns.lineplot(data=data, x=X, y=Y, hue=hue, color=color, marker=marker, markersize=msize, markerfacecolor=mcolor,
                     markeredgecolor="none", linewidth=linewidth, legend=legend, drawstyle=drawstyle, ax=ax)
    else:
        for l in loop:
            if isinstance(l, (int, float)):
                color = colors[l % len(colors)]
            else:
                color = 'lightgray'
            sns.lineplot(data=data[l], x=X, y=Y, hue=hue, color=color, marker=marker, markersize=msize, markerfacecolor=mcolor,
                         markeredgecolor="white", linewidth=linewidth, legend=legend, label=str(l), drawstyle=drawstyle, ax=ax)
    if xlog:
        ax.set(xscale="log")
    if ylog:
        ax.set(yscale="log")

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    if plotname:
        plt.savefig(plotname, bbox_inches='tight', pad_inches=0, transparent=True)
    return ax