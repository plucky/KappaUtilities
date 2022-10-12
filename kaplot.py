# Walter Fontana, 2020

# needs REPAIR

import plotly.graph_objects as go
import matplotlib.pyplot as plt


def show(pdf=''):
    """
    :param pdf: if not '', write to pdf file
    """
    if pdf == '':
        plt.show(block=False)
    else:
        plt.savefig(pdf)


class XY_plot:
    def __init__(self, figsize=(10,8), params={}):
        """
        params: parameter dict to be passed to plot
        """
        self.default_x = 0
        self.default_y = 1
        self.x_axis = None
        self.title = ''
        self.parameters = {'linestyle': '',
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

    def add(self, df, x='', y='', title='', xmajor=0, ymajor=0, params={}):
        """
        ymajor: multiple for major y-tick marks (0 for auto)
        xmajor: multiple for major x-tick marks (0 for auto)
        title: plot title
        params: parameter dict to be passed to plot
        df: pandas dataframe
        y: name of x column
        x: name of y column
        """
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

        arts, = self.ax.plot(df[x], df[y], 'o-', **self.parameters)
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

    snap1 = ks.SnapShot('TestData/snap19.ka')
    sd_df1 = pd.DataFrame(snap1.get_size_distribution(dictionary=True))
    plot = XY_plot(params={'linestyle': '-', 'linewidth': 1., 'markersize': 0})
    plot.add(sd_df1, xmajor=2, ymajor=2000, params={'label': 'snap19', 'color': 'r', 'markerfacecolor': 'r'})
    # show()
    snap2 = ks.SnapShot('TestData/snap98.ka')
    sd_df2 = pd.DataFrame(snap2.get_size_distribution(dictionary=True))
    plot.add(sd_df2, xmajor=2, ymajor=2000, params={'label': 'snap98', 'color': 'g', 'markerfacecolor': 'g'})
    plot.ax.legend()
    show()
