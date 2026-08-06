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


def show(pdf=''):
    """
    :param pdf: if not '', write to pdf file
    """
    if pdf == '':
        plt.show(block=False)
    else:
        plt.savefig(pdf, bbox_inches='tight')


class XYplot:
    def __init__(self, figsize=(10,8), params=None):
        """
        params: parameter dict to be passed to plot
        """
        if not params:
            params = {}
        self.default_x = 0
        self.default_y = 1
        self.x_axis = None
        self.title = ''
        self.parameters = {
                            'linestyle': '',
                            'linewidth': 0.5,
                            'marker': 'o',
                            'label': ''
                           }
        self.parameters = {**self.parameters, **params}
        self.parameters_save = self.parameters
        self.artists = {}
        self.ncurve = 0
        self.legend = None

        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.overlay_axis = None

    def add(self, df, x='', y='', title='', xmajor=0, ymajor=0, params=None):
        """
        ymajor: multiple for major y-tick marks (0 for auto)
        xmajor: multiple for major x-tick marks (0 for auto)
        title: plot title
        params: parameter dict to be passed to plot
        df: pandas dataframe
        y: name of x column
        x: name of y column
        """
        if not params:
            params = {}
        self.parameters = {**self.parameters, **params}
        self.ncurve += 1
        self.parameters['label'] = self.parameters['label'] + f' [{self.ncurve}]'
        if x == '':
            for idx, c in enumerate(df.columns):
                if idx == self.default_x:
                    x = c
                    break
        self.x_axis = x

        if y == '':
            for idx, c in enumerate(df.columns):
                if idx == self.default_y:
                    y = c
                    break
        if title == '':
            self.title = f'{x} vs {y}'
        else:
            self.title = title

        arts, = self.ax.plot(df[x], df[y], **self.parameters)
        self.artists[self.ncurve] = arts

        if xmajor != 0:
            self.ax.xaxis.set_major_locator(plt.MultipleLocator(xmajor))
        if ymajor != 0:
            self.ax.yaxis.set_major_locator(plt.MultipleLocator(ymajor))
        plt.grid(color='lightgrey')
        self.ax.set_xlabel(x)
        self.ax.set_ylabel(y)
        self.ax.set_title(title)
        self.fig.tight_layout()

        self.parameters = self.parameters_save

    def overlay(self, df, y='', ymajor=0, params={}, grid=False):
        """
        grid:
        ymajor: multiple for major y-tick marks (0 for auto)
        params: parameter dict to be passed to plot
        df: pandas dataframe
        y: name of y column
        """
        self.parameters = {**self.parameters, **params}
        self.ncurve += 1
        self.parameters['label'] = self.parameters['label'] + f' [{self.ncurve}]'

        if self.overlay_axis:
            self.overlay_axis.remove()

        self.overlay_axis = self.ax.twinx()  # a second axes that shares the same x-axis

        if y == '':
            for idx, c in enumerate(df.columns):
                if idx == self.default_y:
                    y = c
                    break

        arts, = self.overlay_axis.plot(df[self.x_axis], df[y], 'o-', **self.parameters)
        self.artists[self.ncurve] = arts

        if ymajor != 0:
            self.overlay_axis.yaxis.set_major_locator(plt.MultipleLocator(ymajor))
        if grid:
            plt.grid(color='lightgrey')
        self.overlay_axis.set_ylabel(y)
        self.fig.tight_layout()

        self.parameters = self.parameters_save

    def clear(self, curve_number):
        self.artists[curve_number].remove()
        self.ax.legend()

    def show_legend(self):
        self.legend = self.ax.legend()

    def remove_legend(self):
        if self.legend is not None:
            self.legend.remove()


if __name__ == '__main__':
    import pandas as pd
    import kasnap as ks

    snap1 = ks.SnapShot('TestData/snap__1773.ka')
    sd_df1 = pd.DataFrame(snap1.get_size_distribution(dictionary=True))
    plot = XYplot(params={'linestyle': '-', 'linewidth': 1., 'markersize': 0})
    plot.add(sd_df1, xmajor=20, ymajor=100, params={'label': 'snap1773', 'color': 'r', 'markerfacecolor': 'r'})
    # show()
    snap2 = ks.SnapShot('TestData/snap__1784.ka')
    sd_df2 = pd.DataFrame(snap2.get_size_distribution(dictionary=True))
    plot.add(sd_df2, xmajor=20, ymajor=100, params={'label': 'snap1784', 'color': 'g', 'markerfacecolor': 'g'})

    plot.ax.legend()
    show()