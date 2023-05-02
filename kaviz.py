# Walter Fontana, 2020

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.collections as artcoll
import traceback
import xlsxwriter as xlsx

import kagraph
import kamol


def show():
    plt.show(block=False)


class Canvas:
    def __init__(self, rows=1, cols=1, size=(10, 10), name=None):
        if name is None:
            (filename, line_number, function_name, text) = traceback.extract_stack()[-2]
            name = text[:text.find('=')].strip()
        self.name = name
        self.figure = None
        self.rows = rows
        self.cols = cols

        # axes is a numpy ndarray
        self.figure, self.axes = plt.subplots(nrows=rows, ncols=cols, figsize=size, squeeze=False, frameon=False)
        for i in range(0, self.rows):
            for j in range(0, self.cols):
                self.axes[i, j].axis("off")
                self.axes[i, j].spines["top"].set_visible(False)
                self.axes[i, j].spines["right"].set_visible(False)
                self.axes[i, j].spines["left"].set_visible(False)
                self.axes[i, j].spines["bottom"].set_visible(False)

    def number_of_plot_areas(self):
        return self.rows * self.cols

    def __del__(self):
        if self.figure is not None:
            plt.close(self.figure)
        # print(f'{self.name} deleted')

    def clear(self, i, j):  # assuming i and j start with 1...
        self.axes[i - 1, j - 1].clear()

    def panel2index(self, i, j):  # starting with 1, 1
        return (i - 1) * self.cols + j

    def index2panel(self, n):  # starting with 1
        i = -(- n // self.cols)
        j = n - self.cols * (i - 1)
        return i, j

    def index2axes(self, n, on=False):
        i, j = self.index2panel(n)
        self.axes_visibility(i, j, on=on)
        return self.axes[i - 1, j - 1]

    def panel2axes(self, i, j, on=False):
        self.axes_visibility(i, j, on=on)
        return self.axes[i - 1, j - 1]

    def axes_visibility(self, i, j, on=False):
        i -= 1
        j -= 1
        if on:
            self.axes[i, j].axis("on")
            self.axes[i, j].spines["top"].set_visible(True)
            self.axes[i, j].spines["right"].set_visible(True)
            self.axes[i, j].spines["left"].set_visible(True)
            self.axes[i, j].spines["bottom"].set_visible(True)
        else:
            self.axes[i, j].axis("off")
            self.axes[i, j].spines["top"].set_visible(False)
            self.axes[i, j].spines["right"].set_visible(False)
            self.axes[i, j].spines["left"].set_visible(False)
            self.axes[i, j].spines["bottom"].set_visible(False)


class Renderer:
    def __init__(self, komplex, prog='neato', node_info=False):
        """
        In establishing a Renderer object, a layout of nodes is triggered.
        Subsequently, various display methods can be invoked.

        'prog': neato, dot, twopi, circo, fdp, sfdp, nop, wc, osage, patchwork
        """

        self.ax = None
        self.figure = None

        # generate the networkx representation of komplex
        self.graph = kagraph.KappaGraph(komplex)
        self.nx_graph = self.graph.nxGraph

        self.node_hover_text = {}
        if node_info:
            for type in komplex.composition.keys():
                self.node_hover_text[type] = []
            for node in self.nx_graph.nodes:
                iface = komplex.agents[node]['iface']
                info = f"<b>{node}</b><br>"
                for site in iface.keys():
                    info += f"{site:>10}  ->  state: {iface[site]['state']:<5} bond: {iface[site]['bond']}<br>"
                self.node_hover_text[self.nx_graph.nodes[node]['type']] += [info[:-4]]

        self.nx_palette = ('c', 'r', 'b', 'g', 'm', 'y', 'k', 'w')
        self.html_palette = ('blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'khaki', 'silver')

        self.nx_options = {'font_size': 11,
                           'font_color': 'white',
                           'font_weight': 'bold',
                           'node_size': 400,
                           'labels': {},
                           'edge_color': 'black',
                           'width': 1
                           }

        # assign colors to node types
        self.type_color = {}
        i = 0
        # fill palette index in order of (descending) frequency
        for typ in komplex.composition.keys():
            self.type_color[typ] = i % len(self.nx_palette)
            i += 1

        self.nx_options['node_color'] = []
        for node in self.nx_graph.nodes:
            self.nx_options['node_color'] += [self.nx_palette[self.type_color[self.nx_graph.nodes[node]['type']]]]
        self.legend_colors = [f'{self.nx_palette[self.type_color[n]]}' for n in self.type_color.keys()]

        # layout
        self.positions = nx.nx_agraph.graphviz_layout(self.nx_graph, prog=prog)

    def __del__(self):
        if self.figure is not None:
            plt.close(self.figure)

    def layout(self, prog='neato'):
        self.positions = nx.nx_agraph.graphviz_layout(self.nx_graph, prog=prog)
        self.nx_options['node_color'] = []
        for node in self.nx_graph.nodes:
            self.nx_options['node_color'] += [self.nx_palette[self.type_color[self.nx_graph.nodes[node]['type']]]]
        self.legend_colors = [f'{self.nx_palette[self.type_color[n]]}' for n in self.type_color.keys()]

    def set_palette(self, palette):
        self.nx_palette = palette

    def set_html_palette(self, palette):
        self.html_palette = palette

    def refresh(self, canvas=None, panel=(1, 1), figure_size=(6, 6)):
        if canvas:
            self.ax = canvas.axes[panel[0] - 1, panel[1] - 1]
            self.figure = canvas.figure
            # we clear the panel since we are drawing the whole network
            self.ax.clear()
        else:
            if self.figure:
                plt.close(self.figure)
            self.figure, self.ax = plt.subplots(figsize=figure_size, frameon=False)
            self.ax.axis("off")
            self.ax.spines["top"].set_visible(False)
            self.ax.spines["right"].set_visible(False)
            self.ax.spines["left"].set_visible(False)
            self.ax.spines["bottom"].set_visible(False)
        nx.draw_networkx(self.nx_graph, pos=self.positions, ax=self.ax, **self.nx_options)

    def render(self,
               canvas=None,
               panel=(1, 1),
               labels='short',  # 'type, 'full', 'short';  say 'no' for no labels
               node_size=40,
               font_size=9,
               line_width=1,
               edge_color='gray',
               legend=False,
               title="",
               title_font_size=10,
               title_color='black',
               figure_size=(6, 6)):
        """
        Render a networkx graph with matplotlib.
        """
        self.nx_options['font_size'] = font_size
        self.nx_options['node_size'] = node_size
        # set labels
        self.nx_options['with_labels'] = True
        if labels == 'type':
            self.nx_options['labels'] = {node: self.nx_graph.nodes[node]['type'] for node in self.nx_graph.nodes}
        elif labels == 'short':
            self.nx_options['labels'] = {node: self.nx_graph.nodes[node]['id'] for node in self.nx_graph.nodes}
        elif labels == 'full':
            self.nx_options['labels'] = {node: self.nx_graph.nodes[node]['type'] + self.nx_graph.nodes[node]['id']
                                         for node in self.nx_graph.nodes}
        else:
            self.nx_options['with_labels'] = False

        self.nx_options['edge_color'] = edge_color
        self.nx_options['width'] = line_width  # edge width

        if canvas:
            self.ax = canvas.axes[panel[0] - 1, panel[1] - 1]
            self.figure = canvas.figure
            # we clear the panel since we are drawing the whole network
            self.ax.clear()
        else:
            if self.figure:
                plt.close(self.figure)
            self.figure, self.ax = plt.subplots(figsize=figure_size, frameon=False)
            self.ax.axis("off")
            self.ax.spines["top"].set_visible(False)
            self.ax.spines["right"].set_visible(False)
            self.ax.spines["left"].set_visible(False)
            self.ax.spines["bottom"].set_visible(False)

        nx.draw_networkx(self.nx_graph, pos=self.positions, ax=self.ax, **self.nx_options)

        # the legend
        if legend:
            items = [Line2D([0, 1], [0, 1], color='white', marker='o', markersize=7, markerfacecolor=clr, linewidth=0)
                     for clr in self.legend_colors]
            labels = [f'{node}' for node in self.type_color.keys()]
            self.ax.legend(items, labels)

        if title:
            self.ax.set_title(title, fontsize=title_font_size, color=title_color)

    def color_edge_lists(self, edge_list=[], line_width=1, edge_color='r'):
        # to unify handling, convert to a list of lists (such as coming from a cycle basis)
        if edge_list:
            if not isinstance(edge_list[0], list):
                edge_list = [edge_list]

        self.delete_edge_lists(edge_list=edge_list)
        # draw requested edges in new style
        self.refresh()
        for list_of_edges in edge_list:
            nx.draw_networkx_edges(self.nx_graph, self.positions, ax=self.ax, edgelist=list_of_edges, width=line_width,
                                   edge_color=edge_color)

    def color_node_list(self, node_list=[], color='b', line_width=2):
        node_color = []
        for node in node_list:
            node_color += [self.nx_palette[self.type_color[self.nx_graph.nodes[node]['type']]]]
        self.refresh()
        nx.draw_networkx_nodes(self.nx_graph, nodelist=node_list, pos=self.positions, ax=self.ax,
                               node_size=self.nx_options['node_size'], node_color=node_color,
                               linewidths=line_width, edgecolors=color)

    def delete_edge_lists(self, edge_list=[]):
        # to unify handling, convert to a list of lists (such as coming from a cycle basis)
        if edge_list:
            if not isinstance(edge_list[0], list):
                edge_list = [edge_list]

        untouched_edges = set([frozenset(e) for e in self.nx_graph.edges()])
        for list_of_edges in edge_list:
            untouched_edges = untouched_edges - set([frozenset(e) for e in list_of_edges])
        remaining_edges = [tuple(x) for x in untouched_edges]

        self.refresh()
        self.remove_all_edges()
        # redraw what is left in old style
        nx.draw_networkx_edges(self.nx_graph, self.positions, ax=self.ax, edgelist=remaining_edges, **self.nx_options)

    def delete_node_list(self, node_list=[]):
        untouched_nodes = set(n for n in self.nx_graph.nodes()) - set(n for n in node_list)
        remaining_nodes = [x for x in untouched_nodes]

        self.refresh()
        self.ax.cla()  # clear the whole figure
        node_color = []
        for node in remaining_nodes:
            node_color += [self.nx_palette[self.type_color[self.nx_graph.nodes[node]['type']]]]
        nx.draw_networkx_nodes(self.nx_graph, nodelist=remaining_nodes, pos=self.positions,
                               ax=self.ax, node_size=self.nx_options['node_size'], node_color=node_color)

        # remove the edges incident on the removed nodes
        e_to_delete = []
        for node in node_list:
            e_to_delete += list(self.nx_graph.edges(node))
        edges = set([frozenset(e) for e in self.nx_graph.edges()]) - set([frozenset(e) for e in e_to_delete])
        remaining_edges = [tuple(x) for x in edges]
        nx.draw_networkx_edges(self.nx_graph, self.positions, ax=self.ax, edgelist=remaining_edges, **self.nx_options)

    def remove_all_edges(self):
        for artist in self.ax.get_children():
            if isinstance(artist, artcoll.LineCollection):
                artist.remove()


def show_ranked_complexes(snapshot, canvas, sort='size', cutoff=3, cols=3, rows=1, prog='neato'):
    """
    Display the ranked complexes of a snapshot.

    snapshot: a snapshot object
    sort: 'size' (default) or 'count'
    cols: # of columns of plots
    rows: # of rows of plots
    prog: layout program
    cutoff: size or count cutoff
    """
    ranking = []
    if sort == 'size':
        ranking = sorted(snapshot.complexes, key=lambda x: x.size, reverse=True)
    elif sort == 'count':
        ranking = sorted(snapshot.complexes, key=lambda x: x.count, reverse=True)

    i = 1
    r = []
    for c in ranking[0:cutoff]:
        r.append(Renderer(c, prog=prog))
        r[-1].render(canvas=canvas,
                     panel=canvas.index2panel(i),
                     labels='no',
                     node_size=20,
                     font_size=9,
                     line_width=1,
                     edge_color='gray',
                     title=sort + ' ' + str(c.size))
        i += 1
    plt.show()
    del r


# wrapper
def plot_complex(molecule,
                 labels='no',
                 node_size=40,
                 font_size=9,
                 line_width=1,
                 edge_color='gray',
                 legend=False,
                 title="",
                 title_font_size=10,
                 title_color='black',
                 figure_size=(6, 6)):

    if type(molecule) is str:
        m = kamol.KappaComplex(molecule)
    else:
        m = molecule
    r = Renderer(m)
    r.render(labels=labels,
             node_size=node_size,
             font_size=font_size,
             line_width=line_width,
             edge_color=edge_color,
             legend=legend,
             title=title,
             title_font_size=title_font_size,
             title_color=title_color,
             figure_size=figure_size)
    # show()
    return r


def make_yfile(kappa=None, filename=None):
    """
    makes an input file (xlsx) for the yEd network editor
    """
    nodes = [n for n in kappa.agents]
    node_types = [kappa.agents[n]['info']['type'] for n in kappa.agents]
    workbook = xlsx.Workbook(filename)
    node_sheet = workbook.add_worksheet('Node List')
    edge_sheet = workbook.add_worksheet('Edge List')
    bold = workbook.add_format({'bold': 1})

    # nodes
    node_sheet.write('A1', 'id', bold)
    node_sheet.write('B1', 'node_type', bold)
    for row in range(1, len(nodes)+1):
        node_sheet.write(row, 0, str(nodes[row-1]))
        node_sheet.write(row, 1, node_types[row-1])

    # edges
    source_list = []
    target_list = []
    interaction_list = []
    interaction_labels = []
    for (a1, s1), (a2, s2) in kappa.bonds:
        info = a1 + '@' + s1 + '--' + a2 + '@' + s2
        source_list.append(a1)
        target_list.append(a2)
        sites = info.split('--')
        agent1, _ = kamol.get_identifier(sites[0].split('@')[0])
        agent2, _ = kamol.get_identifier(sites[1].split('@')[0])
        agent_pair = sorted((agent1, agent2))
        txt = f'{agent_pair[0]}-{agent_pair[1]}'
        interaction_list.append(txt)
        interaction_labels.append(info)
    edge_sheet.write('A1', 'source', bold)
    edge_sheet.write('B1', 'target', bold)
    edge_sheet.write('C1', 'interaction', bold)
    edge_sheet.write('D1', 'label', bold)
    for row in range(1, len(source_list)+1):
        edge_sheet.write(row, 0, source_list[row-1])
        edge_sheet.write(row, 1, target_list[row-1])
        edge_sheet.write(row, 2, interaction_labels[row-1])
        edge_sheet.write(row, 3, interaction_list[row - 1])
    workbook.close()


if __name__ == '__main__':
    import kasnap

    # line = open('TestData/bigly.ka', 'r').read()
    # remove newlines that might occur in the file
    # line = re.sub(r'\n+', ' ', line)
    # create a KappaComplex with whatever assignment of node identifiers arises
    # (that's the normalized=False flag).
    line = "A(l[19] r[.] p[2]), A(l[53] r[19] p[42]), A(l[37] r[53] p[45]), A(l[.] r[37] p[29]), P(a1[3] a2[51] a3[" \
           "29] d[.]), A(l[20] r[.] p[3]), A(l[.] r[20] p[27]), P(a1[27] a2[.] a3[.] d[44]), P(a1[.] a2[.] a3[.] d[" \
           "44]), A(l[.] r[.] p[51]), P(a1[14] a2[.] a3[45] d[22]), A(l[24] r[.] p[14]), A(l[.] r[24] p[50]), " \
           "P(a1[50] a2[.] a3[30] d[.]), A(l[13] r[16] p[30]), A(l[.] r[13] p[21]), P(a1[21] a2[35] a3[1] d[9]), " \
           "A(l[.] r[.] p[35]), A(l[.] r[26] p[1]), A(l[26] r[.] p[32]), P(a1[32] a2[.] a3[.] d[9]), A(l[16] r[54] p[" \
           ".]), A(l[54] r[38] p[40]), A(l[38] r[.] p[7]), P(a1[.] a2[.] a3[7] d[52]), P(a1[18] a2[.] a3[4] d[52]), " \
           "A(l[.] r[.] p[18]), A(l[23] r[.] p[4]), A(l[49] r[23] p[47]), A(l[28] r[49] p[11]), A(l[6] r[28] p[31]), " \
           "A(l[.] r[6] p[12]), P(a1[.] a2[12] a3[.] d[.]), P(a1[.] a2[31] a3[.] d[15]), P(a1[.] a2[5] a3[.] d[15]), " \
           "A(l[.] r[.] p[5]), P(a1[8] a2[.] a3[11] d[36]), A(l[.] r[39] p[8]), A(l[39] r[.] p[25]), P(a1[41] a2[.] " \
           "a3[25] d[33]), A(l[.] r[.] p[41]), P(a1[.] a2[43] a3[.] d[33]), A(l[.] r[46] p[43]), A(l[46] r[.] p[10]), " \
           "P(a1[10] a2[.] a3[.] d[22]), P(a1[34] a2[.] a3[.] d[36]), A(l[.] r[.] p[34]), P(a1[.] a2[47] a3[.] d[.]), " \
           "P(a1[.] a2[40] a3[42] d[17]), P(a1[48] a2[2] a3[.] d[17]), A(l[.] r[.] p[48]) "
    line2 = "A(l[.] r[4] p[1]), A(l[4] r[.] p[3]), P(a1[3] a2[1] a3[.] d[2]), P(a1[.] a2[.] a3[.] d[2])"

    # argument can be a KappaMolecule or a string
    r = plot_complex(line)
    del r

    c1 = kamol.KappaComplex(line)
    r1 = Renderer(c1)
    r1.render(labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')

    print(c1)
    print(c1.canonical)

    c2 = kamol.Canonical2Expression(c1.canonical, c1.system_views)
    r2 = Renderer(c2)
    r2.render(labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')

    plt.ion()
    canvas = Canvas(2, 1)
    r1.render(canvas, panel=(1, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # input()
    r2.render(canvas, panel=(2, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')

    plt.ion()
    canvas = Canvas(3, 1)
    snap = kasnap.SnapShot('TestData/snap__1773.ka')
    show_ranked_complexes(snap, canvas=canvas)

    kappa_ring = 'A(r[6] l[1]),A(r[1] l[2]),A(r[2] l[3]),A(r[3] l[4]),A(r[4] l[5]),A(r[5] l[6])'
    c = kamol.KappaComplex(kappa_ring)
    print(c.show())
    r = Renderer(c)
    r.render(node_size=200)
    g = kagraph.KappaGraph(c)
    cycle = g.get_cycle()
    print(cycle)

    # r.new_plot()
    # r.color_edge_lists(edge_list=[cycle[:-1]], line_width=5, edge_color='red')
    # show()
    # r.delete_edge_lists(edge_list=[cycle])
    # show()
    # r.render(node_size=200)
    # show()
