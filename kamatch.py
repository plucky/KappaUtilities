# Walter Fontana, 2020-2022
"""
Graph embedding by exploiting site graph rigidity; efficient and straightforward.
"""

from collections import deque


# The graphs are given by the internal representation in 'kamol.KappaMolecule',
# which, for the purpose of this code, behaves like a networkx graph.

def print_map(maps):
    for i in range(0, len(maps)):
        print(f'map {i + 1}:')
        for k, v in maps[i].items():
            print(f'{k} --> {v}')

# ================================================================================================
# The rigid way... (site graphs only)
# ================================================================================================


class Fail(Exception):
    pass


# In the following implementation, speed-ups were achieved by
# * avoiding assignment of host and pattern graphs (not sure why that has an effect)
# * turning off degree checking as preflight
# * avoiding creating a new SiteGraphMatcher object at each comparison in all_embeddings
#   (hence all matching functions are now methods of the matcher object)


class SiteGraphMatcher:
    """
    Implements graph pattern matching for Kappa exploiting the rigidity of site graphs.
    Does handle multi-graphs, but not graphs with multiple disconnected components.
    """

    def __init__(self):
        self.p_start = ''
        self.h_start = ''
        # tentative subgraph embedding; if pattern lacks a bond that host has, it's OK
        self.sub = False
        self.mapping = {}

    def isomorphic(self, host, pattern):
        # Check size
        if host.size != pattern.size:
            return False
        # Check composition
        if host.sum_formula != pattern.sum_formula:
            return False
        # Check local properties; turned off because of overhead
        # d1 = sorted(d for n, d in host.degree())
        # d2 = sorted(d for n, d in pattern.degree())
        # if d1 != d2:
        #     return False
        return self.morphism(host, pattern)

    def automorphisms(self, pattern):
        # wrapper for special case of all_embeddings()
        # start, stop = pattern.type_slice[pattern.rarest_type]
        mappings = []
        # for node in pattern.name_list[start:stop]:
        for lst in pattern.type_slice:
            for node in lst:
                if self._embed(pattern, pattern, h_start=node):
                    # sort sensibly for readability
                    self.mapping = {k: v for k, v in sorted(self.mapping.items(), key=lambda x: x[0])}
                    # the rigidity approach never produces identical embeddings
                    mappings += [self.mapping]
        return mappings

    def isomorphic_unsafe(self, host, pattern, h_start, p_start):
        # Here we supply specific anchors for the match, rather than checking all possible
        # anchors as in (safe) isomorphic().

        # Obviously, check size first
        if host.size != pattern.size:
            return False
        # Check composition
        if host.sum_formula != pattern.sum_formula:
            return False
        # Check local properties; turned off because of overhead
        # d1 = sorted(d for n, d in host.degree())
        # d2 = sorted(d for n, d in pattern.degree())
        # if d1 != d2:
        #     return False
        return self._embed(host, pattern, h_start=h_start, p_start=p_start)

    def embed(self, host, pattern):
        # Check size
        if host.size < pattern.size:
            return False
        # Check composition
        for node_type in pattern.composition:
            if node_type not in host.composition.keys():
                return False
            if host.composition[node_type] < pattern.composition[node_type]:
                return False
        return self.morphism(host, pattern)

    def sub_embed(self, host, pattern):
        self.sub = True
        answer = self.embed(host, pattern)
        self.sub = False
        return answer

    def morphism(self, host, pattern):
        # start, stop = host.type_slice[pattern.rarest_type]
        # for node in host.name_list[start:stop]:
        for lst in host.type_slice:
            for node in lst:
                if self._embed(host, pattern, h_start=node):
                    return True
        return False

    def all_embeddings(self, host, pattern):
        # Check size
        if host.size < pattern.size:
            return []
        # Check composition
        for node_type in pattern.composition:
            if node_type not in host.composition.keys():
                return []
            if host.composition[node_type] < pattern.composition[node_type]:
                return []

        # start, stop = host.type_slice[pattern.rarest_type]
        mappings = []
        # for node in host.name_list[start:stop]:
        for lst in host.type_slice:
            for node in lst:
                if self._embed(host, pattern, h_start=node):
                    # sort sensibly for readability
                    self.mapping = {k: v for k, v in sorted(self.mapping.items(), key=lambda x: x[0])}
                    # the rigidity approach never produces identical embeddings
                    mappings += [self.mapping]
        return mappings

    def number_of_all_embeddings(self, host, pattern):
        # Check size
        if host.size < pattern.size:
            return 0
        # Check composition
        for node_type in pattern.composition:
            if node_type not in host.composition.keys():
                return 0
            if host.composition[node_type] < pattern.composition[node_type]:
                return 0
        num = 0
        for lst in host.type_slice:
            for node in lst:
                if self._embed(host, pattern, h_start=node):
                    num += 1
        return num

    def _embed(self, host, pattern, p_start=None, h_start=None):

        if not p_start:
            self.p_start = pattern.embedding_anchor
        else:
            self.p_start = p_start
        if not h_start:
            self.h_start = host.embedding_anchor
        else:
            self.h_start = h_start
        self.mapping = {}

        try:
            self.traverse(host, pattern, self.p_start, self.h_start)
            return True
        except Fail:
            return False

    def traverse(self, host, pattern, p_node, h_node):
        # sans recursion
        stack = deque()
        visited = set()
        # stack (current pattern node, predecessor, current host node)
        stack.append((p_node, '', h_node))
        while stack:
            p_node, parent_p_node, h_node = stack.pop()
            if p_node not in visited:
                if parent_p_node:
                    # the site at which we left the last pattern node to reach the current pattern node
                    site = pattern.navigation[(parent_p_node, p_node)]
                    # the agent that is bound on that site but in the host graph
                    # (it doesn't matter if there are multiple bonds between two nodes; all we care about
                    # is reaching the node in host that must be matched against the node in pattern.)
                    h_node, x, x = host.agents[h_node]['iface'][site]['bond'].partition(pattern.bond_sep)
                if self.node_match(host, pattern, h_node, p_node):
                    # update the mapping
                    self.mapping[p_node] = h_node
                    visited.add(p_node)
                    for neighbor in pattern[p_node]:
                        # stack the neighbors
                        stack.append((neighbor, p_node, h_node))
                else:
                    raise Fail

    # def traverse_with_recursion(self, p_node, h_node):
    #     self.p_visited[p_node] = True
    #     if self.start:
    #         self.start = False
    #     else:
    #         # the site at which we left the last pattern node to reach the current pattern node
    #         last_p_node = list(self.stack)[-1]  # peek
    #         site = self.pattern.navigation[(last_p_node, p_node)]
    #         # the agent that is bound on that site but in the host graph
    #         h_node = self.host.agents[h_node]['iface'][site]['bond'].split(self.pattern.bond_sep)[0]
    #     if not self.node_match(h_node, p_node):
    #         raise Fail
    #     else:
    #         # update the mapping
    #         self.mapping[p_node] = h_node
    #         # store the last p_node
    #         self.stack.append(p_node)
    #         for neighbor in self.pattern[p_node]:
    #             if not self.p_visited[neighbor]:
    #                 self.traverse_with_recursion(neighbor, h_node)
    #         self.stack.pop()

    def node_match_bonds_only(self, host, pattern, h_node, p_node):
        # no state comparison for binding-only models !
        #
        # the prefix 'h' stands for 'host' and 'p' for 'pattern'
        # start with type match
        h_node_type = host.agents[h_node]['info']['type']
        p_node_type = pattern.agents[p_node]['info']['type']
        if h_node_type != p_node_type:
            return False

        h_node_iface = host.agents[h_node]['iface']
        p_node_iface = pattern.agents[p_node]['iface']

        for site_name in p_node_iface:
            if site_name not in h_node_iface:
                return False
            else:
                h_bond = h_node_iface[site_name]['bond']
                p_bond = p_node_iface[site_name]['bond']

                if p_bond == '.':
                    if h_bond != '.':
                        return False
                elif p_bond == '#':  # don't care
                    return True
                elif p_bond == '_':  # unspecific bond
                    if h_bond == '.':
                        return False
                else:  # specific bond
                    if h_bond == '.':
                        return False
                    else:
                        # both sites are bound
                        p_partner, _, p_site = p_bond.partition('@')
                        h_partner, _, h_site = h_bond.partition('@')
                        if p_partner in self.mapping:
                            if not (self.mapping[p_partner] == h_partner):
                                return False
                        if h_site != p_site:
                            return False
        return True

    # general version
    def node_match(self, host, pattern, h_node, p_node):
        # here, the prefix 'h' stands for 'host' and 'p' for 'pattern'
        # start with type match
        h_node_type = host.agents[h_node]['info']['type']
        p_node_type = pattern.agents[p_node]['info']['type']
        if h_node_type != p_node_type:
            return False

        h_node_iface = host.agents[h_node]['iface']
        p_node_iface = pattern.agents[p_node]['iface']

        for site_name in p_node_iface:
            if site_name not in h_node_iface:
                return False
            else:
                if p_node_iface[site_name]['state'] != '#':
                    if p_node_iface[site_name]['state'] != h_node_iface[site_name]['state']:
                        return False

                h_bond = h_node_iface[site_name]['bond']
                p_bond = p_node_iface[site_name]['bond']
                if '@' in h_bond:
                    h_partner, x, h_site = h_bond.partition(host.bond_sep)
                if '@' in p_bond:
                    p_partner, x, p_site = p_bond.partition(pattern.bond_sep)

                if p_bond == '.':
                    if h_bond != '.':
                        if not self.sub:
                            return False
                elif '@' in p_bond:  # specific bond
                    if not ('@' in h_bond):
                        return False
                    else:
                        # both sites are bound
                        if p_partner in self.mapping:
                            if not (self.mapping[p_partner] == h_partner):
                                return False
                        if h_site != p_site:
                            return False
                elif p_bond == '_':  # unspecific bond
                    if h_bond == '.':
                        return False
                elif '.' in p_bond:  # stub ('.', as in free, is caught above)
                    if h_bond == '.':  # the site is free
                        return False
                    elif h_bond == '_':
                        return False  # the pattern state is more specific (stub) than the host state
                    elif '@' in h_bond:
                        ghost_site, ghost_type = p_bond.split('.')
                        h_type = h_partner.split('.')[0]
                        if ghost_type != h_type or ghost_site != h_site:
                            return False
                    elif '.' in h_bond:  # h_bond is also a stub
                        if p_bond != h_bond:
                            return False
        return True

# -------------------------------------------------------------------------------------------


if __name__ == '__main__':
    import kamol
    import time
    import re

    # SGM = SiteGraphMatcher()
    # G1 = kamol.kappaObject('B(b[1]{a}), B(b[1])')
    # print(G1.show())
    # maps = SGM.automorphisms(G1)
    # print(f'{len(maps)} automorphisms:')
    # print_map(maps)
    # print(f'{"-"*20}')
    # G1 = kamol.kappaObject('A(a[1],b[2]),A(a[1],b[3]),B(a[2],c[4]),B(a[3],c[5]),C(b[4],x{u}), C(b[5],x{p})')
    # G1.show()
    # maps = SGM.automorphisms(G1)
    # print(f'{len(maps)} automorphisms:')
    # print_map(maps)

    SGM = SiteGraphMatcher()
    G1 = kamol.KappaComplex('A(a[1],b[2]),A(a[1],b[3]),B(a[2],c[4]),B(a[3],c[5]),C(b[4],x{u}[.]),C(b[5],x{p}[.])')
    G2 = kamol.KappaComplex('A(a[1],b[2]),A(a[1],b[3]),B(a[2],c[4]),B(a[3],c[5]),C(b[5],x{p}[.]),C(b[4],x{u}[.])')
    print(G1.show())
    print(G2.show())
    maps = SGM.automorphisms(G1)
    print(f'{len(maps)} automorphisms:')
    print_map(maps)
    maps = SGM.automorphisms(G2)
    print(f'{len(maps)} automorphisms:')
    print_map(maps)
    print()
    print(G1.canonical)
    print(G2.canonical)

    # SGM = SiteGraphMatcher()
    # G1 = kamol.kappaObject('B(b[1]), B(b[1]{a})')
    # G2 = kamol.kappaObject('B(b[2]{b}), B(b[2]{a})')
    # maps = SGM.all_embeddings(G1, G2)
    # print_map(maps)
    # maps = SGM.all_embeddings(G2, G1)
    # print_map(maps)
    #
    # # G1 = kamol.kappaObject('A(r[5] p[.]), A(l[1] p[2]), A(r[1] p[3] l[5]), P(a3[2] d[4]), P(a1[3] d[4])')
    # # G2a = kamol.kappaObject('A(l[1] p[2]), A(r[1] p[3]), P(a3[2] d[4]), P(a1[3] d[4])')
    # # G2c = kamol.kappaObject('A(l[1] p[2]), A(r[1] p[3]), P(a1[2] d[4]), P(a1[3] d[4])')
    # # G2b = kamol.kappaObject('A(l[1] p[2]), A(r[1] p[3]), P(a2[2] d[4]), P(a3[3] d[4])')
    #
    # G1 = kamol.kappaObject('A(b[1] a[2]), A(b[3] a[2]), B(a[1] x{p}), B(a[3] x{u})')
    # G1.show()
    # G2 = kamol.kappaObject('A(b[3] a[2]), A(b[1] a[2]), B(a[1] x{p}), B(a[3] x{u})')
    # G2.show()
    #
    # SGM = SiteGraphMatcher()
    # print(f'rigid: G1 and G2 are isomorphic: {SGM.isomorphic_unsafe(G1, G2, None, None)}')
    # print(f'rigid: G1 and G2 are isomorphic: {SGM.isomorphic_unsafe(G1, G2, "A.1.", "A.2.")}')
    # print(f'rigid: G1 and G2 are isomorphic: {SGM.isomorphic(G1, G2)}')
    #
    # # usage scenarios
    #
    # G1 = kamol.kappaObject('A(b[1] a[2]), A(b[3] a[2]), B(a[1] x{p}), B(a[3] x{u})')
    # G2 = kamol.kappaObject('A(b[2]), B(a[2] x{u})')
    # G1.show()
    # G2.show()
    # SGM = SiteGraphMatcher()
    # maps = SGM.all_embeddings(G1, G2)
    # print(f'number of embeddings of G2 into G1: {len(maps)} ')
    # print_map(maps)
    #
    # # input a file containing one (large) Kappa string
    # line = open('TestData/bigly.ka', 'r').read()
    # # remove newlines that might occur in the file
    # line = re.sub(r'\n+', ' ', line)
    # # create a KappaComplex with whatever assignment of node identifiers arises
    # # (that's the normalize=False flag).
    # G1 = kamol.kappaObject(line)
    # G2 = kamol.kappaObject(line)
    # G2.show()
    # print(f'G1 and G2 are isomorphic: {SGM.isomorphic(G1, G2)}')
    #
    # G1 = kamol.kappaObject('A(x[1],y[2]), B(x[2],y[3]), C(x[3],y[1])')
    # G2 = kamol.kappaObject('A(x[1],y[2]), B(x[2],y[3]), C(x[3],y[1])')
    # for i in range(0, 10):
    #     # randomly permute identifiers:
    #     G2.remap(change='randomize')
    #     print(f'G1 and G2 are isomorphic: {SGM.isomorphic(G1, G2)}')
    # print('')
    #
    # G1 = kamol.kappaObject('A(x[1],y[2]), B(x[2],y[3]), C(x[3],y[1])')
    # G2 = kamol.kappaObject('A(x[.],y[2]), B(x[2],y[3]), C(x[3],y[.])')
    # print('G1:')
    # G1.show()
    # print('G2:')
    # G2.show()
    # print(f'G1 and G2 are isomorphic: {SGM.isomorphic(G1, G2)}')
    # print(f'G2 is embeddable in G1: {SGM.isomorphic(G1, G2)}')
    #
    # G1 = kamol.kappaObject('A(x[1],y[2]), A(x[1],y[3]), C(x[2]), C(x[3]{p})')
    # G2 = kamol.kappaObject('A(x[1],y[2]), A(x[1],y[3]), C(x[2]), C(x[3])')
    # print('G1:')
    # G1.show()
    # print('G2:')
    # G2.show()
    # print(f'G2 is embeddable in G1: {SGM.all_embeddings(G1, G2)}')
    # print(f'G1 and G2 are isomorphic: {SGM.all_embeddings(G1, G2)}')
    #
    # G1 = kamol.kappaObject('A(x[.] y[2]), A(x[2] y[3]), A(x[3] y[4]), A(x[4] y[1]), B(x[1])')
    # G2 = kamol.kappaObject('A(x y[2]), A(x[2] y)')
    # print('G1:')
    # G1.show()
    # print('G2:')
    # G2.show()
    # maps = SGM.all_embeddings(G1, G2)
    # print(f'number of embeddings from G2 into G1: {len(maps)} ')
    # print_map(maps)
    #
    # G1 = kamol.kappaObject('A(b[1] a[2]), A(b[3] a[2]), B(a[1] x{p}), B(a[3] x{u})')
    # G2 = kamol.kappaObject('A(b[2]), B(a[2] x)')
    # maps = SGM.all_embeddings(G1, G2)
    # print(f'number of embeddings of G2 into G1: {len(maps)} ')
    # print_map(maps)
    # #
    # G1 = kamol.kappaObject('A(b[1] a[2]), A(b[1] a[2], c[3]), B(a[3] x[.]{p})')
    # G2 = kamol.kappaObject('A(b[1]), A(b[1])')
    # G1.show()
    # G2.show()
    # maps = SGM.all_embeddings(G1, G2)
    # print(f'number of embeddings of G2 into G1: {len(maps)} ')
    # print_map(maps)
    #
    # G1 = kamol.kappaObject('A(b[1] a[2]), A(b[3] a[2]), B(a[1] x{p}), B(a[3] x{u})')
    # G2 = kamol.kappaObject('A(a[1]), A(a[1])')
    # G1.show()
    # G2.show()
    # maps = SGM.all_embeddings(G1, G2)
    # print(f'number of embeddings of G2 into G1: {len(maps)} ')
    # print_map(maps)
    # #
    # G1 = kamol.kappaObject('A(x[1] y[2]), A(x[2] y[3]), A(x[3] y[4]]), A(x[4] y[1]})')
    # G2 = kamol.kappaObject('A(x[1] y[2]), A(x[2] y[3]), A(x[3] y[4]]), A(x[4] y[1]})')
    # maps = SGM.all_embeddings(G1, G2)
    # print(f'number of embeddings of G2 into G1: {len(maps)}')
    # print_map(maps)
    #
    # # --------------------------------------------------------------------------------
    #
    # print("big stuff:")
    # line = open('TestData/bigly.ka', 'r').read()
    # line = re.sub(r'\n+', ' ', line)
    # start = time.process_time()
    # c1 = kamol.kappaObject(line)
    # end = time.process_time()
    # print(f'seconds: {end - start}')
    # start = time.process_time()
    # print(SGM.isomorphic(c1, c1))
    # end = time.process_time()
    # print(f'seconds: {end - start}')
