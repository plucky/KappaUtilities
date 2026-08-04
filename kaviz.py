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

    def color_edge_lists(self, edge_list=None, line_width=1, edge_color='r'):
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

    def color_node_list(self, node_list=None, color='b', line_width=2):
        node_color = []
        for node in node_list:
            node_color += [self.nx_palette[self.type_color[self.nx_graph.nodes[node]['type']]]]
        self.refresh()
        nx.draw_networkx_nodes(self.nx_graph, nodelist=node_list, pos=self.positions, ax=self.ax,
                               node_size=self.nx_options['node_size'], node_color=node_color,
                               linewidths=line_width, edgecolors=color)

    def delete_edge_lists(self, edge_list=None):
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

    def delete_node_list(self, node_list=None):
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
        m = kamol.kappa_to_molecule(molecule)
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
    # line = "A(l[19] r[.] p[2]), A(l[53] r[19] p[42]), A(l[37] r[53] p[45]), A(l[.] r[37] p[29]), P(a1[3] a2[51] a3[" \
    #        "29] d[.]), A(l[20] r[.] p[3]), A(l[.] r[20] p[27]), P(a1[27] a2[.] a3[.] d[44]), P(a1[.] a2[.] a3[.] d[" \
    #        "44]), A(l[.] r[.] p[51]), P(a1[14] a2[.] a3[45] d[22]), A(l[24] r[.] p[14]), A(l[.] r[24] p[50]), " \
    #        "P(a1[50] a2[.] a3[30] d[.]), A(l[13] r[16] p[30]), A(l[.] r[13] p[21]), P(a1[21] a2[35] a3[1] d[9]), " \
    #        "A(l[.] r[.] p[35]), A(l[.] r[26] p[1]), A(l[26] r[.] p[32]), P(a1[32] a2[.] a3[.] d[9]), A(l[16] r[54] p[" \
    #        ".]), A(l[54] r[38] p[40]), A(l[38] r[.] p[7]), P(a1[.] a2[.] a3[7] d[52]), P(a1[18] a2[.] a3[4] d[52]), " \
    #        "A(l[.] r[.] p[18]), A(l[23] r[.] p[4]), A(l[49] r[23] p[47]), A(l[28] r[49] p[11]), A(l[6] r[28] p[31]), " \
    #        "A(l[.] r[6] p[12]), P(a1[.] a2[12] a3[.] d[.]), P(a1[.] a2[31] a3[.] d[15]), P(a1[.] a2[5] a3[.] d[15]), " \
    #        "A(l[.] r[.] p[5]), P(a1[8] a2[.] a3[11] d[36]), A(l[.] r[39] p[8]), A(l[39] r[.] p[25]), P(a1[41] a2[.] " \
    #        "a3[25] d[33]), A(l[.] r[.] p[41]), P(a1[.] a2[43] a3[.] d[33]), A(l[.] r[46] p[43]), A(l[46] r[.] p[10]), " \
    #        "P(a1[10] a2[.] a3[.] d[22]), P(a1[34] a2[.] a3[.] d[36]), A(l[.] r[.] p[34]), P(a1[.] a2[47] a3[.] d[.]), " \
    #        "P(a1[.] a2[40] a3[42] d[17]), P(a1[48] a2[2] a3[.] d[17]), A(l[.] r[.] p[48]) "
    # line2 = "A(l[.] r[4] p[1]), A(l[4] r[.] p[3]), P(a1[3] a2[1] a3[.] d[2]), P(a1[.] a2[.] a3[.] d[2])"
    #
    # # argument can be a KappaMolecule or a string
    # r = plot_complex(line)
    # del r

    line = "B(d[1] a1[.] a2[2] a3[3]), A(l[4] r[5] b[6]), A(l[7] r[8] b[9]), B(d[10] a1[.] a2[11] a3[12]), B(d[13] a1[14] a2[.] a3[15]), A(l[16] r[17] b[18]), A(l[19] r[20] b[21]), A(l[22] r[.] b[23]), A(l[24] r[25] b[26]), A(l[27] r[28] b[29]), A(l[30] r[.] b[31]), B(d[.] a1[32] a2[33] a3[.]), A(l[34] r[35] b[36]), A(l[37] r[38] b[39]), A(l[.] r[40] b[41]), B(d[42] a1[.] a2[43] a3[44]), A(l[45] r[46] b[47]), A(l[48] r[49] b[50]), A(l[51] r[52] b[53]), B(d[54] a1[.] a2[.] a3[47]), A(l[.] r[55] b[56]), A(l[57] r[58] b[59]), B(d[60] a1[61] a2[62] a3[63]), A(l[64] r[65] b[66]), A(l[67] r[68] b[69]), A(l[70] r[71] b[72]), B(d[73] a1[74] a2[.] a3[.]), A(l[75] r[76] b[77]), B(d[78] a1[79] a2[80] a3[81]), A(l[82] r[83] b[84]), B(d[85] a1[86] a2[.] a3[87]), A(l[88] r[.] b[89]), A(l[90] r[57] b[91]), A(l[92] r[93] b[94]), B(d[95] a1[96] a2[29] a3[.]), A(l[97] r[98] b[99]), B(d[100] a1[.] a2[101] a3[.]), B(d[102] a1[103] a2[104] a3[41]), B(d[105] a1[.] a2[106] a3[107]), A(l[108] r[109] b[110]), B(d[111] a1[112] a2[113] a3[114]), B(d[115] a1[.] a2[116] a3[117]), A(l[118] r[119] b[120]), A(l[121] r[122] b[123]), A(l[.] r[124] b[125]), B(d[126] a1[127] a2[128] a3[129]), A(l[130] r[131] b[132]), A(l[.] r[133] b[134]), B(d[135] a1[136] a2[.] a3[137]), B(d[138] a1[139] a2[140] a3[141]), A(l[8] r[142] b[143]), B(d[144] a1[.] a2[145] a3[.]), A(l[146] r[147] b[148]), B(d[149] a1[150] a2[.] a3[151]), A(l[65] r[152] b[153]), A(l[.] r[154] b[155]), B(d[156] a1[.] a2[157] a3[69]), A(l[158] r[159] b[160]), A(l[161] r[162] b[163]), A(l[164] r[165] b[166]), A(l[167] r[168] b[169]), A(l[170] r[171] b[172]), A(l[173] r[174] b[175]), A(l[176] r[177] b[178]), A(l[179] r[130] b[3]), B(d[.] a1[180] a2[91] a3[181]), A(l[182] r[183] b[184]), B(d[185] a1[.] a2[.] a3[186]), A(l[187] r[.] b[188]), A(l[189] r[190] b[191]), A(l[192] r[193] b[194]), A(l[162] r[195] b[196]), B(d[197] a1[198] a2[.] a3[199]), A(l[200] r[201] b[202]), A(l[203] r[204] b[205]), A(l[206] r[207] b[208]), B(d[209] a1[210] a2[211] a3[212]), A(l[.] r[213] b[214]), B(d[215] a1[216] a2[.] a3[217]), B(d[218] a1[219] a2[220] a3[221]), A(l[222] r[.] b[223]), B(d[.] a1[224] a2[225] a3[.]), B(d[226] a1[227] a2[228] a3[229]), B(d[230] a1[231] a2[.] a3[.]), A(l[232] r[233] b[234]), A(l[235] r[236] b[237]), A(l[.] r[238] b[239]), A(l[124] r[.] b[240]), A(l[25] r[241] b[242]), B(d[243] a1[175] a2[244] a3[245]), B(d[246] a1[247] a2[31] a3[.]), A(l[248] r[.] b[249]), B(d[250] a1[251] a2[252] a3[253]), A(l[254] r[255] b[256]), A(l[257] r[258] b[259]), A(l[260] r[261] b[262]), A(l[263] r[264] b[265]), A(l[.] r[266] b[267]), B(d[268] a1[9] a2[.] a3[269]), A(l[270] r[271] b[272]), A(l[273] r[274] b[275]), A(l[276] r[82] b[277]), B(d[149] a1[278] a2[.] a3[.]), A(l[279] r[7] b[280]), A(l[281] r[282] b[283]), A(l[284] r[285] b[286]), A(l[171] r[164] b[287]), B(d[288] a1[.] a2[.] a3[289]), B(d[290] a1[.] a2[.] a3[39]), A(l[291] r[292] b[293]), B(d[294] a1[295] a2[296] a3[297]), A(l[298] r[299] b[300]), A(l[301] r[302] b[303]), A(l[304] r[305] b[306]), B(d[307] a1[308] a2[.] a3[309]), A(l[310] r[311] b[312]), A(l[313] r[314] b[315]), A(l[316] r[317] b[318]), B(d[319] a1[320] a2[.] a3[.]), A(l[321] r[322] b[323]), B(d[324] a1[.] a2[325] a3[326]), A(l[327] r[328] b[329]), A(l[330] r[331] b[332]), B(d[333] a1[334] a2[335] a3[336]), A(l[337] r[338] b[339]), B(d[226] a1[340] a2[341] a3[.]), B(d[342] a1[.] a2[343] a3[344]), A(l[345] r[346] b[347]), B(d[348] a1[349] a2[.] a3[.]), A(l[350] r[351] b[352]), B(d[353] a1[354] a2[.] a3[.]), A(l[355] r[356] b[357]), B(d[.] a1[89] a2[358] a3[359]), A(l[282] r[360] b[361]), A(l[362] r[363] b[364]), B(d[365] a1[.] a2[366] a3[249]), A(l[367] r[368] b[369]), A(l[370] r[371] b[372]), A(l[373] r[51] b[374]), A(l[.] r[375] b[376]), A(l[377] r[235] b[378]), B(d[379] a1[380] a2[.] a3[.]), A(l[.] r[381] b[228]), A(l[382] r[383] b[384]), B(d[385] a1[.] a2[386] a3[.]), A(l[387] r[19] b[79]), B(d[388] a1[.] a2[389] a3[.]), A(l[390] r[391] b[392]), B(d[393] a1[394] a2[.] a3[.]), A(l[119] r[395] b[396]), B(d[397] a1[398] a2[399] a3[400]), B(d[401] a1[.] a2[.] a3[.]), A(l[.] r[45] b[62]), B(d[402] a1[.] a2[403] a3[.]), B(d[404] a1[405] a2[.] a3[406]), A(l[407] r[408] b[114]), A(l[409] r[410] b[411]), B(d[412] a1[.] a2[.] a3[413]), A(l[414] r[415] b[416]), A(l[417] r[418] b[419]), A(l[420] r[421] b[422]), A(l[423] r[424] b[425]), B(d[426] a1[.] a2[427] a3[428]), A(l[429] r[430] b[431]), B(d[.] a1[.] a2[.] a3[432]), A(l[433] r[.] b[434]), B(d[435] a1[.] a2[.] a3[436]), A(l[58] r[420] b[437]), B(d[438] a1[.] a2[439] a3[440]), A(l[441] r[442] b[443]), B(d[54] a1[.] a2[.] a3[.]), A(l[444] r[445] b[446]), A(l[447] r[448] b[449]), B(d[450] a1[451] a2[.] a3[452]), A(l[453] r[454] b[405]), B(d[455] a1[456] a2[.] a3[457]), A(l[458] r[30] b[459]), A(l[460] r[461] b[462]), B(d[463] a1[.] a2[464] a3[.]), A(l[415] r[465] b[466]), A(l[467] r[67] b[468]), A(l[469] r[470] b[471]), B(d[156] a1[.] a2[472] a3[473]), A(l[474] r[475] b[476]), A(l[363] r[.] b[477]), B(d[478] a1[479] a2[480] a3[481]), B(d[412] a1[482] a2[483] a3[484]), A(l[485] r[.] b[486]), A(l[487] r[488] b[199]), A(l[133] r[489] b[490]), B(d[455] a1[491] a2[492] a3[188]), A(l[49] r[493] b[296]), A(l[494] r[355] b[413]), B(d[.] a1[495] a2[496] a3[497]), A(l[498] r[499] b[500]), B(d[501] a1[267] a2[.] a3[502]), A(l[503] r[504] b[74]), B(d[505] a1[506] a2[507] a3[508]), A(l[509] r[510] b[511]), B(d[512] a1[312] a2[26] a3[.]), A(l[513] r[514] b[515]), B(d[516] a1[.] a2[517] a3[.]), A(l[152] r[518] b[519]), A(l[520] r[182] b[150]), B(d[450] a1[521] a2[522] a3[.]), A(l[523] r[257] b[136]), A(l[524] r[525] b[220]), A(l[526] r[284] b[406]), A(l[527] r[528] b[529]), B(d[530] a1[.] a2[531] a3[.]), A(l[532] r[533] b[534]), A(l[535] r[536] b[537]), B(d[102] a1[538] a2[539] a3[540]), B(d[294] a1[425] a2[541] a3[542]), A(l[543] r[544] b[545]), A(l[546] r[547] b[548]), B(d[549] a1[.] a2[550] a3[223]), B(d[551] a1[552] a2[.] a3[490]), A(l[.] r[553] b[554]), B(d[1] a1[555] a2[556] a3[557]), A(l[558] r[559] b[560]), A(l[.] r[561] b[562]), A(l[563] r[564] b[565]), A(l[566] r[567] b[568]), B(d[569] a1[570] a2[300] a3[571]), A(l[317] r[572] b[573]), A(l[190] r[498] b[574]), B(d[575] a1[.] a2[132] a3[.]), B(d[576] a1[.] a2[577] a3[.]), A(l[17] r[578] b[579]), B(d[580] a1[581] a2[315] a3[582]), A(l[583] r[108] b[584]), B(d[246] a1[293] a2[.] a3[585]), A(l[.] r[586] b[587]), A(l[588] r[589] b[590]), A(l[591] r[592] b[593]), A(l[46] r[232] b[594]), A(l[98] r[595] b[596]), B(d[.] a1[597] a2[598] a3[599]), A(l[600] r[601] b[602]), A(l[204] r[603] b[604]), B(d[605] a1[.] a2[606] a3[607]), B(d[397] a1[202] a2[608] a3[609]), A(l[610] r[611] b[612]), A(l[613] r[614] b[542]), A(l[615] r[616] b[617]), B(d[618] a1[619] a2[.] a3[.]), B(d[.] a1[620] a2[172] a3[372]), A(l[.] r[621] b[432]), A(l[622] r[623] b[624]), B(d[625] a1[626] a2[.] a3[627]), B(d[628] a1[629] a2[.] a3[596]), A(l[630] r[631] b[632]), A(l[633] r[634] b[635]), A(l[636] r[637] b[607]), B(d[638] a1[639] a2[640] a3[.]), A(l[641] r[642] b[643]), A(l[644] r[.] b[225]), B(d[645] a1[646] a2[647] a3[648]), A(l[649] r[650] b[651]), A(l[652] r[653] b[556]), B(d[654] a1[655] a2[.] a3[656]), A(l[299] r[657] b[658]), B(d[218] a1[659] a2[660] a3[.]), A(l[71] r[661] b[662]), A(l[76] r[663] b[427]), B(d[365] a1[664] a2[665] a3[.]), A(l[666] r[667] b[627]), B(d[668] a1[560] a2[669] a3[670]), A(l[653] r[671] b[457]), B(d[353] a1[155] a2[672] a3[.]), A(l[572] r[.] b[673]), A(l[674] r[675] b[676]), A(l[174] r[.] b[677]), B(d[678] a1[679] a2[680] a3[681]), B(d[682] a1[.] a2[683] a3[612]), A(l[684] r[310] b[685]), A(l[686] r[687] b[688]), B(d[689] a1[.] a2[.] a3[.]), A(l[690] r[691] b[692]), A(l[693] r[16] b[106]), A(l[694] r[695] b[672]), A(l[696] r[447] b[452]), B(d[697] a1[256] a2[529] a3[698]), A(l[699] r[.] b[107]), B(d[700] a1[701] a2[632] a3[702]), A(l[703] r[704] b[701]), A(l[.] r[390] b[521]), A(l[705] r[370] b[670]), B(d[706] a1[707] a2[708] a3[.]), A(l[709] r[666] b[710]), A(l[711] r[712] b[713]), A(l[331] r[699] b[714]), A(l[715] r[64] b[716]), B(d[105] a1[717] a2[718] a3[272]), B(d[719] a1[720] a2[431] a3[721]), B(d[385] a1[.] a2[.] a3[.]), A(l[475] r[722] b[112]), B(d[250] a1[723] a2[724] a3[725]), A(l[726] r[.] b[11]), B(d[727] a1[.] a2[728] a3[729]), B(d[730] a1[731] a2[732] a3[733]), A(l[734] r[.] b[735]), B(d[736] a1[134] a2[.] a3[737]), A(l[.] r[738] b[739]), A(l[740] r[260] b[741]), B(d[.] a1[.] a2[.] a3[742]), A(l[743] r[744] b[745]), B(d[746] a1[747] a2[.] a3[.]), A(l[748] r[749] b[217]), B(d[750] a1[.] a2[.] a3[751]), A(l[20] r[752] b[753]), A(l[754] r[755] b[723]), A(l[756] r[686] b[757]), A(l[758] r[759] b[760]), A(l[761] r[703] b[762]), A(l[763] r[764] b[765]), A(l[766] r[767] b[552]), B(d[750] a1[768] a2[769] a3[770]), A(l[771] r[772] b[773]), B(d[774] a1[775] a2[.] a3[194]), B(d[776] a1[242] a2[.] a3[153]), A(l[777] r[778] b[779]), B(d[780] a1[781] a2[.] a3[.]), A(l[255] r[782] b[783]), B(d[426] a1[593] a2[.] a3[.]), A(l[784] r[785] b[786]), A(l[787] r[526] b[788]), A(l[789] r[.] b[309]), B(d[790] a1[791] a2[287] a3[745]), A(l[792] r[793] b[794]), B(d[795] a1[.] a2[.] a3[796]), A(l[.] r[600] b[797]), A(l[798] r[.] b[648]), B(d[799] a1[800] a2[.] a3[.]), A(l[.] r[801] b[802]), A(l[803] r[22] b[804]), A(l[.] r[805] b[366]), B(d[806] a1[807] a2[.] a3[.]), A(l[808] r[809] b[810]), B(d[697] a1[811] a2[812] a3[651]), A(l[504] r[813] b[814]), A(l[815] r[158] b[816]), B(d[115] a1[.] a2[817] a3[.]), A(l[818] r[819] b[820]), A(l[236] r[821] b[394]), B(d[516] a1[.] a2[.] a3[275]), A(l[601] r[822] b[724]), B(d[746] a1[762] a2[166] a3[823]), B(d[824] a1[.] a2[802] a3[825]), A(l[826] r[827] b[828]), A(l[829] r[830] b[831]), A(l[832] r[535] b[833]), B(d[834] a1[602] a2[.] a3[835]), B(d[836] a1[56] a2[.] a3[837]), B(d[838] a1[.] a2[.] a3[239]), A(l[518] r[566] b[839]), A(l[442] r[674] b[721]), B(d[138] a1[.] a2[.] a3[.]), B(d[268] a1[840] a2[841] a3[.]), A(l[842] r[377] b[843]), A(l[844] r[845] b[731]), B(d[846] a1[847] a2[.] a3[848]), A(l[849] r[330] b[531]), A(l[755] r[850] b[851]), A(l[852] r[853] b[854]), B(d[230] a1[.] a2[855] a3[.]), A(l[856] r[313] b[857]), B(d[625] a1[858] a2[859] a3[.]), A(l[860] r[861] b[473]), B(d[333] a1[862] a2[.] a3[259]), A(l[863] r[429] b[113]), A(l[864] r[161] b[791]), B(d[85] a1[865] a2[866] a3[867]), A(l[391] r[868] b[869]), A(l[870] r[871] b[480]), A(l[872] r[873] b[874]), B(d[875] a1[277] a2[.] a3[876]), A(l[877] r[387] b[878]), B(d[879] a1[214] a2[880] a3[.]), A(l[445] r[881] b[127]), A(l[882] r[815] b[835]), B(d[126] a1[883] a2[.] a3[486]), A(l[884] r[885] b[664]), B(d[886] a1[.] a2[.] a3[887]), B(d[.] a1[888] a2[.] a3[889]), B(d[890] a1[110] a2[891] a3[892]), A(l[893] r[726] b[894]), A(l[895] r[.] b[812]), B(d[618] a1[896] a2[.] a3[.]), A(l[897] r[898] b[899]), A(l[900] r[.] b[901]), A(l[.] r[787] b[902]), B(d[719] a1[.] a2[.] a3[.]), A(l[165] r[633] b[464]), B(d[903] a1[904] a2[635] a3[905]), B(d[906] a1[833] a2[907] a3[.]), A(l[908] r[417] b[451]), B(d[909] a1[.] a2[.] a3[.]), A(l[564] r[.] b[910]), B(d[911] a1[.] a2[912] a3[.]), A(l[913] r[914] b[915]), A(l[448] r[362] b[43]), B(d[73] a1[.] a2[878] a3[352]), A(l[916] r[917] b[647]), A(l[.] r[918] b[919]), B(d[920] a1[921] a2[922] a3[.]), A(l[923] r[884] b[539]), A(l[924] r[925] b[491]), B(d[926] a1[927] a2[928] a3[929]), B(d[930] a1[931] a2[466] a3[932]), A(l[933] r[.] b[934]), B(d[935] a1[.] a2[.] a3[.]), A(l[661] r[936] b[61]), B(d[185] a1[.] a2[.] a3[.]), A(l[691] r[453] b[181]), A(l[.] r[.] b[656]), B(d[903] a1[937] a2[938] a3[765]), A(l[578] r[939] b[940]), A(l[809] r[281] b[941]), B(d[942] a1[.] a2[.] a3[943]), B(d[944] a1[945] a2[946] a3[.]), A(l[738] r[844] b[86]), A(l[947] r[948] b[889]), A(l[949] r[950] b[951]), A(l[952] r[953] b[728]), A(l[954] r[955] b[956]), B(d[957] a1[.] a2[396] a3[958]), A(l[959] r[960] b[961]), A(l[962] r[963] b[964]), B(d[965] a1[966] a2[378] a3[.]), A(l[967] r[.] b[968]), B(d[790] a1[163] a2[.] a3[810]), A(l[147] r[969] b[970]), A(l[603] r[829] b[971]), A(l[.] r[734] b[581]), A(l[663] r[972] b[973]), B(d[974] a1[975] a2[.] a3[976]), A(l[873] r[.] b[977]), A(l[614] r[978] b[63]), A(l[979] r[316] b[245]), A(l[980] r[981] b[349]), A(l[266] r[613] b[768]), B(d[580] a1[.] a2[.] a3[982]), A(l[.] r[24] b[983]), B(d[935] a1[.] a2[.] a3[984]), A(l[985] r[986] b[198]), A(l[987] r[988] b[989]), B(d[774] a1[990] a2[991] a3[.]), A(l[.] r[842] b[992]), B(d[342] a1[662] a2[123] a3[.]), A(l[168] r[949] b[993]), A(l[994] r[696] b[608]), B(d[926] a1[.] a2[84] a3[995]), B(d[996] a1[449] a2[.] a3[.]), B(d[.] a1[997] a2[998] a3[476]), B(d[999] a1[1000] a2[.] a3[1001]), A(l[1002] r[1003] b[577]), B(d[1004] a1[854] a2[.] a3[459]), A(l[1005] r[1006] b[117]), B(d[605] a1[1007] a2[.] a3[548]), A(l[493] r[1008] b[1009]), A(l[1010] r[532] b[1011]), A(l[1012] r[222] b[483]), B(d[.] a1[1013] a2[.] a3[973]), B(d[1014] a1[.] a2[1015] a3[.]), A(l[744] r[118] b[1016]), A(l[1017] r[1018] b[1019]), A(l[1020] r[1021] b[180]), A(l[383] r[1022] b[598]), A(l[1023] r[771] b[247]), A(l[969] r[444] b[140]), B(d[824] a1[.] a2[604] a3[1024]), A(l[657] r[.] b[655]), A(l[1025] r[1026] b[1027]), A(l[.] r[1028] b[912]), A(l[1029] r[1030] b[1031]), A(l[1032] r[527] b[937]), A(l[1033] r[985] b[210]), A(l[.] r[48] b[1034]), B(d[.] a1[.] a2[1035] a3[1036]), A(l[.] r[1037] b[1038]), A(l[1039] r[1040] b[269]), A(l[898] r[872] b[87]), B(d[875] a1[.] a2[.] a3[468]), B(d[1041] a1[191] a2[.] a3[.]), A(l[.] r[740] b[1042]), A(l[936] r[206] b[1043]), B(d[1044] a1[.] a2[820] a3[.]), A(l[1045] r[1046] b[186]), A(l[.] r[1047] b[380]), A(l[1048] r[1049] b[1050]), B(d[628] a1[1051] a2[.] a3[.]), A(l[621] r[.] b[1052]), A(l[1053] r[1054] b[1055]), B(d[957] a1[283] a2[562] a3[1056]), B(d[974] a1[1057] a2[.] a3[.]), A(l[1058] r[1059] b[729]), A(l[1060] r[690] b[932]), B(d[1061] a1[.] a2[1062] a3[.]), A(l[1063] r[.] b[931]), B(d[1064] a1[.] a2[1065] a3[1066]), A(l[1067] r[487] b[128]), B(d[1068] a1[1069] a2[.] a3[.]), A(l[1070] r[1071] b[211]), A(l[1072] r[200] b[1073]), B(d[1074] a1[.] a2[237] a3[437]), A(l[1075] r[1076] b[1077]), A(l[675] r[947] b[1078]), A(l[258] r[1079] b[1080]), B(d[1044] a1[53] a2[.] a3[.]), B(d[307] a1[739] a2[983] a3[.]), A(l[764] r[588] b[456]), B(d[654] a1[1081] a2[.] a3[.]), A(l[948] r[1082] b[1083]), A(l[1084] r[1085] b[946]), B(d[1086] a1[.] a2[919] a3[1087]), A(l[1088] r[1045] b[1089]), B(d[463] a1[.] a2[.] a3[.]), A(l[1090] r[1039] b[1091]), A(l[360] r[.] b[1092]), A(l[1093] r[1094] b[896]), B(d[393] a1[.] a2[303] a3[940]), A(l[1095] r[1012] b[640]), B(d[197] a1[.] a2[.] a3[.]), B(d[1096] a1[.] a2[753] a3[1097]), A(l[650] r[959] b[1098]), A(l[.] r[1099] b[867]), A(l[1100] r[1101] b[141]), A(l[667] r[1102] b[1103]), B(d[1104] a1[332] a2[.] a3[1105]), A(l[1106] r[276] b[775]), B(d[1107] a1[1108] a2[1109] a3[.]), A(l[586] r[273] b[1110]), B(d[.] a1[329] a2[.] a3[.]), A(l[430] r[863] b[497]), A(l[533] r[1111] b[585]), A(l[.] r[485] b[400]), A(l[1071] r[179] b[12]), B(d[549] a1[.] a2[794] a3[.]), A(l[559] r[1112] b[104]), B(d[319] a1[.] a2[831] a3[1113]), A(l[1114] r[849] b[1115]), B(d[1116] a1[.] a2[1117] a3[416]), A(l[885] r[.] b[1118]), A(l[1119] r[967] b[1120]), A(l[1121] r[.] b[538]), B(d[1122] a1[.] a2[1123] a3[.]), A(l[822] r[1124] b[796]), B(d[1125] a1[961] a2[50] a3[.]), A(l[845] r[1126] b[145]), B(d[1127] a1[.] a2[773] a3[1128]), A(l[375] r[994] b[1129]), A(l[950] r[1130] b[1131]), B(d[575] a1[.] a2[.] a3[1132]), A(l[1133] r[1032] b[680]), A(l[418] r[1134] b[847]), B(d[1135] a1[369] a2[6] a3[.]), A(l[960] r[1048] b[837]), A(l[233] r[962] b[876]), B(d[838] a1[1136] a2[.] a3[1042]), A(l[1137] r[877] b[1128]), B(d[438] a1[.] a2[446] a3[.]), A(l[1138] r[1088] b[609]), A(l[1139] r[1140] b[597]), A(l[801] r[1141] b[1065]), A(l[408] r[1142] b[1143]), B(d[1086] a1[617] a2[1144] a3[.]), A(l[1102] r[248] b[1097]), A(l[381] r[1145] b[1146]), A(l[1147] r[1148] b[683]), B(d[512] a1[169] a2[1019] a3[1149]), A(l[1150] r[1151] b[1152]), B(d[700] a1[.] a2[519] a3[23]), A(l[.] r[.] b[1153]), B(d[215] a1[1154] a2[1155] a3[1027]), A(l[1156] r[1157] b[229]), A(l[1158] r[520] b[995]), B(d[638] a1[.] a2[234] a3[36]), A(l[28] r[1159] b[278]), B(d[1160] a1[.] a2[658] a3[910]), A(l[1161] r[1162] b[399]), A(l[528] r[1163] b[506]), A(l[813] r[1020] b[1164]), B(d[1165] a1[1166] a2[72] a3[1167]), B(d[1168] a1[1169] a2[843] a3[1170]), B(d[846] a1[1171] a2[.] a3[.]), A(l[623] r[298] b[101]), A(l[410] r[895] b[1172]), B(d[1173] a1[851] a2[1143] a3[1103]), A(l[939] r[1060] b[1001]), B(d[1165] a1[1115] a2[.] a3[1131]), A(l[917] r[433] b[1174]), B(d[1175] a1[1080] a2[1176] a3[.]), A(l[1177] r[1178] b[905]), B(d[1168] a1[1179] a2[.] a3[.]), A(l[1180] r[1114] b[1181]), A(l[544] r[1182] b[1183]), A(l[1112] r[458] b[702]), A(l[963] r[954] b[340]), A(l[778] r[1017] b[927]), A(l[1054] r[1184] b[770]), A(l[371] r[1093] b[96]), B(d[1185] a1[.] a2[1186] a3[.]), A(l[305] r[1187] b[1066]), A(l[1188] r[1158] b[938]), A(l[122] r[754] b[571]), B(d[930] a1[1189] a2[1190] a3[537]), A(l[.] r[1023] b[646]), A(l[368] r[1191] b[1192]), A(l[.] r[.] b[1166]), A(l[547] r[1193] b[343]), A(l[1194] r[1029] b[103]), A(l[1148] r[1195] b[1196]), A(l[1140] r[.] b[1176]), A(l[1191] r[367] b[698]), B(d[1197] a1[1198] a2[915] a3[1199]), B(d[834] a1[1181] a2[1200] a3[99]), A(l[356] r[1090] b[1189]), A(l[925] r[1201] b[811]), A(l[292] r[1202] b[1203]), A(l[1204] r[709] b[224]), A(l[1193] r[1205] b[389]), A(l[1206] r[748] b[769]), A(l[1207] r[1208] b[619]), A(l[759] r[641] b[732]), B(d[1209] a1[1210] a2[783] a3[1172]), A(l[561] r[1161] b[1211]), B(d[42] a1[1212] a2[673] a3[.]), A(l[1213] r[254] b[44]), A(l[1214] r[1215] b[665]), B(d[95] a1[.] a2[.] a3[941]), A(l[1184] r[1216] b[507]), A(l[.] r[321] b[479]), B(d[111] a1[376] a2[94] a3[120]), A(l[1217] r[1218] b[921]), A(l[782] r[.] b[1219]), B(d[909] a1[902] a2[.] a3[21]), B(d[886] a1[1220] a2[786] a3[.]), A(l[819] r[615] b[582]), A(l[595] r[1221] b[129]), B(d[1209] a1[1073] a2[574] a3[18]), A(l[1222] r[1223] b[1224]), B(d[1225] a1[1226] a2[1227] a3[.]), A(l[1228] r[1121] b[1229]), B(d[60] a1[1230] a2[.] a3[.]), A(l[1231] r[291] b[1232]), A(l[1233] r[1234] b[550]), A(l[1126] r[1235] b[1186]), B(d[478] a1[.] a2[735] a3[.]), B(d[505] a1[1236] a2[814] a3[1237]), A(l[.] r[1238] b[1239]), A(l[1208] r[.] b[33]), B(d[1240] a1[1011] a2[.] a3[.]), A(l[1241] r[1242] b[892]), A(l[642] r[705] b[1243]), B(d[1244] a1[.] a2[1243] a3[1245]), A(l[1246] r[1247] b[639]), A(l[1248] r[1058] b[1249]), B(d[1240] a1[1250] a2[.] a3[1251]), B(d[1096] a1[1083] a2[.] a3[1252]), A(l[1182] r[1253] b[1254]), B(d[1255] a1[1256] a2[964] a3[857]), A(l[.] r[423] b[439]), A(l[1257] r[1258] b[599]), A(l[1259] r[1177] b[998]), B(d[1260] a1[1261] a2[.] a3[1262]), A(l[1263] r[1067] b[1264]), A(l[351] r[1084] b[289]), B(d[100] a1[992] a2[951] a3[1265]), A(l[1266] r[509] b[880]), A(l[553] r[636] b[1267]), A(l[1268] r[.] b[733]), B(d[1269] a1[.] a2[306] a3[1270]), A(l[1238] r[1214] b[555]), A(l[1163] r[1271] b[717]), A(l[955] r[1272] b[825]), A(l[1130] r[1266] b[253]), A(l[1273] r[1274] b[334]), A(l[868] r[.] b[1275]), B(d[1276] a1[374] a2[.] a3[1038]), A(l[1277] r[558] b[482]), A(l[1278] r[1279] b[1280]), A(l[274] r[1281] b[484]), B(d[668] a1[1203] a2[.] a3[1282]), A(l[1218] r[270] b[1283]), A(l[611] r[1284] b[1081]), A(l[201] r[1285] b[325]), A(l[.] r[583] b[1286]), A(l[986] r[1072] b[707]), B(d[402] a1[.] a2[.] a3[.]), B(d[290] a1[1232] a2[.] a3[677]), A(l[1287] r[897] b[336]), A(l[1288] r[337] b[945]), A(l[1289] r[27] b[344]), B(d[78] a1[1290] a2[1291] a3[1292]), A(l[1284] r[893] b[1256]), A(l[499] r[543] b[1198]), A(l[1281] r[1293] b[1294]), A(l[159] r[1295] b[32]), B(d[1296] a1[.] a2[1297] a3[.]), A(l[853] r[167] b[1132]), B(d[1298] a1[.] a2[779] a3[901]), B(d[1175] a1[.] a2[1211] a3[1299]), A(l[1300] r[1301] b[1302]), A(l[1006] r[1303] b[1051]), A(l[1304] r[818] b[1305]), A(l[.] r[173] b[1306]), A(l[93] r[1188] b[1123]), A(l[1124] r[.] b[502]), A(l[5] r[1307] b[1220]), A(l[461] r[860] b[1308]), A(l[207] r[187] b[251]), A(l[1195] r[1075] b[495]), A(l[1234] r[146] b[1265]), A(l[1059] r[1150] b[883]), A(l[1253] r[1309] b[997]), A(l[1310] r[.] b[14]), A(l[1311] r[170] b[841]), B(d[920] a1[.] a2[.] a3[160]), B(d[1312] a1[993] a2[692] a3[1313]), A(l[510] r[1314] b[1315]), A(l[.] r[513] b[1236]), A(l[1316] r[1317] b[1318]), A(l[634] r[467] b[1319]), A(l[1295] r[1228] b[212]), A(l[972] r[1095] b[718]), B(d[1244] a1[.] a2[.] a3[.]), A(l[1317] r[1320] b[481]), A(l[1008] r[382] b[865]), B(d[1260] a1[.] a2[804] a3[.]), B(d[135] a1[590] a2[.] a3[1091]), A(l[1321] r[1287] b[326]), B(d[551] a1[.] a2[1092] a3[419]), A(l[421] r[1322] b[1108]), A(l[.] r[882] b[1117]), A(l[827] r[1222] b[116]), B(d[996] a1[1323] a2[.] a3[.]), A(l[314] r[630] b[626]), A(l[871] r[263] b[80]), A(l[1003] r[761] b[1324]), B(d[1122] a1[1325] a2[1239] a3[573]), A(l[767] r[1326] b[888]), A(l[.] r[693] b[81]), A(l[981] r[524] b[1105]), B(d[1014] a1[.] a2[545] a3[.]), A(l[1028] r[1248] b[1327]), A(l[40] r[1328] b[1290]), A(l[1046] r[1231] b[1329]), B(d[1330] a1[.] a2[1077] a3[1331]), A(l[1111] r[1332] b[1170]), B(d[1160] a1[1324] a2[579] a3[816]), B(d[1333] a1[280] a2[.] a3[1334]), A(l[687] r[591] b[227]), B(d[1312] a1[.] a2[.] a3[788]), A(l[1307] r[1304] b[428]), B(d[799] a1[.] a2[.] a3[1335]), A(l[1301] r[1336] b[1337]), A(l[1099] r[1338] b[1323]), A(l[1339] r[1340] b[522]), B(d[806] a1[.] a2[1341] a3[1280]), A(l[1322] r[.] b[1270]), A(l[988] r[1321] b[1292]), A(l[695] r[908] b[1342]), A(l[1343] r[1133] b[1056]), B(d[1104] a1[240] a2[565] a3[.]), A(l[1344] r[441] b[297]), A(l[1258] r[1100] b[1345]), A(l[1346] r[373] b[1347]), B(d[906] a1[.] a2[1348] a3[.]), A(l[1293] r[1349] b[606]), B(d[379] a1[.] a2[1034] a3[956]), A(l[.] r[.] b[751]), A(l[1350] r[1351] b[1262]), A(l[1352] r[327] b[252]), A(l[1353] r[1063] b[1354]), B(d[569] a1[1355] a2[.] a3[1055]), A(l[850] r[789] b[1245]), A(l[.] r[1138] b[1356]), B(d[1074] a1[515] a2[.] a3[1016]), A(l[.] r[1357] b[660]), A(l[322] r[1358] b[1000]), B(d[879] a1[1359] a2[.] a3[148]), A(l[285] r[1217] b[1149]), B(d[795] a1[1315] a2[.] a3[.]), A(l[1360] r[92] b[976]), B(d[999] a1[1098] a2[1152] a3[1183]), A(l[1361] r[1362] b[216]), A(l[1358] r[1363] b[1036]), A(l[1351] r[407] b[320]), B(d[1127] a1[125] a2[.] a3[1305]), A(l[38] r[649] b[1226]), A(l[1040] r[1364] b[1365]), B(d[144] a1[1366] a2[.] a3[.]), A(l[712] r[1310] b[1169]), B(d[388] a1[716] a2[1367] a3[.]), A(l[1368] r[1311] b[800]), A(l[1142] r[1300] b[1369]), B(d[1370] a1[59] a2[.] a3[760]), A(l[142] r[1025] b[231]), B(d[965] a1[1371] a2[1043] a3[1365]), A(l[589] r[1268] b[943]), A(l[1274] r[176] b[1372]), A(l[793] r[1373] b[1057]), A(l[261] r[1002] b[472]), A(l[1364] r[1374] b[966]), A(l[978] r[763] b[866]), B(d[944] a1[1089] a2[77] a3[970]), A(l[1134] r[1375] b[1210]), A(l[536] r[.] b[517]), A(l[1030] r[952] b[1035]), A(l[1223] r[777] b[496]), B(d[1068] a1[1376] a2[.] a3[.]), B(d[1225] a1[1129] a2[1377] a3[713]), A(l[1378] r[715] b[817]), B(d[324] a1[.] a2[624] a3[205]), A(l[1326] r[301] b[308]), B(d[706] a1[797] a2[1379] a3[.]), A(l[1340] r[.] b[840]), A(l[183] r[1194] b[1167]), B(d[736] a1[.] a2[1153] a3[968]), A(l[1380] r[1273] b[492]), A(l[1381] r[1278] b[1382]), A(l[1141] r[1380] b[922]), A(l[35] r[1263] b[341]), A(l[830] r[684] b[929]), B(d[1333] a1[.] a2[.] a3[1383]), A(l[271] r[121] b[887]), A(l[821] r[923] b[15]), A(l[264] r[4] b[1251]), A(l[1384] r[1241] b[1250]), B(d[1064] a1[.] a2[1385] a3[511]), A(l[1386] r[1339] b[1199]), A(l[525] r[1070] b[139]), B(d[836] a1[.] a2[1283] a3[1196]), B(d[1185] a1[.] a2[422] a3[934]), A(l[722] r[1005] b[1387]), A(l[195] r[694] b[1366]), A(l[1388] r[1206] b[1212]), A(l[616] r[980] b[540]), A(l[1037] r[192] b[984]), B(d[1389] a1[977] a2[1390] a3[1391]), A(l[1082] r[1033] b[1390]), B(d[645] a1[710] a2[1392] a3[1393]), A(l[514] r[523] b[221]), B(d[530] a1[584] a2[1394] a3[434]), A(l[.] r[1204] b[1395]), A(l[302] r[856] b[508]), A(l[1336] r[758] b[737]), A(l[1328] r[.] b[1261]), A(l[338] r[784] b[958]), A(l[346] r[345] b[1154]), A(l[1349] r[.] b[1007]), B(d[1255] a1[462] a2[.] a3[.]), A(l[.] r[1257] b[137]), A(l[861] r[34] b[848]), B(d[682] a1[477] a2[.] a3[411]), A(l[1162] r[743] b[398]), B(d[1389] a1[568] a2[.] a3[.]), A(l[1094] r[987] b[386]), A(l[1396] r[711] b[1325]), A(l[637] r[1397] b[742]), B(d[209] a1[1395] a2[714] a3[1318]), A(l[1235] r[.] b[907]), A(l[109] r[1259] b[659]), A(l[1338] r[1119] b[1190]), A(l[1320] r[652] b[858]), A(l[1279] r[1398] b[1399]), B(d[1330] a1[.] a2[208] a3[676]), A(l[1047] r[1400] b[1109]), B(d[1269] a1[1302] a2[1050] a3[1372]), A(l[1247] r[1139] b[1355]), A(l[1272] r[1316] b[1348]), A(l[1285] r[.] b[295]), B(d[348] a1[.] a2[.] a3[1401]), A(l[1402] r[1207] b[1403]), A(l[1357] r[832] b[1404]), A(l[1159] r[756] b[1405]), B(d[1135] a1[1406] a2[1387] a3[685]), B(d[243] a1[1249] a2[.] a3[66]), A(l[1397] r[798] b[620]), B(d[1407] a1[1219] a2[384] a3[184]), A(l[1408] r[97] b[1144]), A(l[1363] r[1381] b[1409]), A(l[55] r[.] b[1252]), B(d[942] a1[.] a2[262] a3[688]), A(l[1026] r[644] b[1391]), B(d[1107] a1[.] a2[286] a3[1337]), A(l[1314] r[.] b[1410]), A(l[1411] r[409] b[1376]), A(l[1332] r[.] b[747]), A(l[395] r[1010] b[725]), A(l[1271] r[1343] b[1230]), A(l[1201] r[1412] b[1394]), A(l[.] r[1053] b[2]), A(l[1398] r[1289] b[1413]), B(d[890] a1[1329] a2[1414] a3[1118]), A(l[1415] r[1396] b[244]), A(l[1309] r[90] b[403]), A(l[328] r[913] b[1416]), A(l[1375] r[1156] b[1383]), B(d[911] a1[1405] a2[1404] a3[.]), A(l[1417] r[350] b[904]), A(l[567] r[916] b[1299]), A(l[1187] r[70] b[335]), A(l[1418] r[622] b[1334]), B(d[1419] a1[.] a2[.] a3[.]), A(l[914] r[1420] b[681]), B(d[1173] a1[339] a2[1421] a3[443]), A(l[1157] r[1277] b[807]), A(l[1422] r[88] b[1237]), A(l[465] r[503] b[1313]), B(d[404] a1[.] a2[1120] a3[.]), A(l[1412] r[.] b[1379]), A(l[1423] r[864] b[862]), A(l[1362] r[1415] b[440]), A(l[1303] r[1213] b[1406]), A(l[881] r[1417] b[1069]), A(l[918] r[546] b[982]), A(l[631] r[75] b[1392]), A(l[1151] r[1423] b[157]), A(l[704] r[826] b[1200]), A(l[671] r[1346] b[891]), B(d[1061] a1[196] a2[.] a3[.]), A(l[1221] r[1422] b[436]), A(l[1216] r[1288] b[1385]), B(d[1407] a1[318] a2[.] a3[894]), A(l[785] r[1402] b[629]), A(l[1178] r[1353] b[1421]), A(l[1021] r[1233] b[541]), A(l[1101] r[766] b[359]), A(l[1018] r[563] b[1015]), A(l[752] r[189] b[1335]), B(d[776] a1[757] a2[392] a3[1354]), A(l[68] r[414] b[1087]), A(l[454] r[1344] b[975]), B(d[1116] a1[1052] a2[143] a3[1275]), B(d[435] a1[.] a2[1413] a3[869]), A(l[489] r[1352] b[1136]), A(l[.] r[1246] b[1424]), A(l[1145] r[1360] b[781]), A(l[592] r[1386] b[1359]), B(d[727] a1[1294] a2[.] a3[989]), A(l[1374] r[.] b[1291]), A(l[.] r[304] b[859]), A(l[470] r[469] b[1155]), B(d[401] a1[265] a2[1399] a3[1425]), B(d[10] a1[1009] a2[.] a3[.]), A(l[.] r[1361] b[219]), A(l[1426] r[870] b[669]), A(l[131] r[852] b[1013]), B(d[1041] a1[594] a2[361] a3[534]), A(l[1373] r[1427] b[557]), A(l[1242] r[203] b[990]), B(d[730] a1[1254] a2[1409] a3[1267]), B(d[678] a1[1229] a2[1306] a3[1403]), A(l[772] r[1147] b[1341]), A(l[1420] r[1368] b[1377]), A(l[213] r[1378] b[570]), A(l[177] r[924] b[1401]), A(l[1085] r[.] b[1171]), A(l[.] r[1384] b[358]), A(l[311] r[979] b[823]), B(d[1419] a1[741] a2[.] a3[471]), A(l[1205] r[1137] b[991]), A(l[.] r[1411] b[1297]), B(d[1370] a1[1382] a2[364] a3[1345]), A(l[1215] r[1180] b[1179]), A(l[424] r[494] b[1393]), A(l[83] r[1106] b[855]), A(l[1076] r[1426] b[1113]), B(d[1125] a1[1192] a2[347] a3[874]), A(l[805] r[279] b[1367]), A(l[238] r[1350] b[928]), B(d[13] a1[.] a2[1410] a3[.]), B(d[.] a1[1110] a2[643] a3[.]), A(l[52] r[1418] b[1371]), B(d[1276] a1[1327] a2[839] a3[1078]), A(l[953] r[.] b[1062]), B(d[689] a1[.] a2[.] a3[1319]), A(l[193] r[37] b[1414]), B(d[1428] a1[971] a2[899] a3[500]), A(l[241] r[1388] b[151]), B(d[780] a1[1224] a2[1342] a3[1146]), A(l[1427] r[933] b[1331]), B(d[501] a1[.] a2[1347] a3[1424]), B(d[1197] a1[1429] a2[.] a3[554]), A(l[749] r[803] b[720]), B(d[1296] a1[.] a2[1031] a3[1164]), A(l[154] r[808] b[708]), B(d[1428] a1[1369] a2[178] a3[.]), B(d[576] a1[.] a2[323] a3[1174]), A(l[1400] r[.] b[1425]), A(l[1202] r[474] b[1227]), A(l[488] r[1408] b[1024]), B(d[288] a1[357] a2[1308] a3[.]), A(l[.] r[610] b[354]), A(l[1049] r[460] b[1282]), B(d[1298] a1[1264] a2[828] a3[587]), A(l[1022] r[792] b[1429]), B(d[1004] a1[1356] a2[1416] a3[1286]), A(l[1079] r[900] b[679])"
    c1 = kamol.kappa_to_molecule(line)
    r1 = Renderer(c1)
    r1.render(labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')

    # make_yfile(kappa=c1, filename='/Users/wf7/Desktop/Temp/test.xlsx')
    # print(c1)
    # print(c1.canonical)
    #
    # c2 = kamol.canonical_to_representation(c1.canonical, c1.local_view_index)
    # r2 = Renderer(c2)
    # r2.render(labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    #
    # plt.ion()
    # canvas = Canvas(2, 1)
    # r1.render(canvas, panel=(1, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # # input()
    # r2.render(canvas, panel=(2, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')

    # plt.ion()
    # canvas = Canvas(3, 1)
    # snap = kasnap.SnapShot('TestData/snap__1773.ka')
    # show_ranked_complexes(snap, canvas=canvas)
    #
    # kappa_ring = 'A(r[6] l[1]),A(r[1] l[2]),A(r[2] l[3]),A(r[3] l[4]),A(r[4] l[5]),A(r[5] l[6])'
    # c = kamol.kappa_to_molecule(kappa_ring)
    # print(c.show())
    # r = Renderer(c)
    # r.render(node_size=200)
    # g = kagraph.KappaGraph(c)
    # cycle = g.get_cycle()
    # print(cycle)

    # r.new_plot()
    # r.color_edge_lists(edge_list=[cycle[:-1]], line_width=5, edge_color='red')
    # show()
    # r.delete_edge_lists(edge_list=[cycle])
    # show()
    # r.render(node_size=200)
    # show()
