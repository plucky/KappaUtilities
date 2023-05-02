# Walter Fontana, 2022
"""
This module defines the internal representation of a kappa molecule. It implements the parser, and a
 number of utilities, such as copying molecules, binding two molecules, removing a bond from a molecule,
 and identifying fragmentation.
"""
import re
import sys
import random
import pprint
import ujson
from collections import deque
from collections import defaultdict


class ParseFail(Exception):
    pass


def convert(text): return int(text) if text.isdigit() else text


def alphanum_key(key): return [convert(c) for c in re.split('([0-9]+)', key)]


def bond2type(b):
    x = sorted([(re.sub(r'.\d+.', '', b[0][0]), b[0][1]), (re.sub(r'.\d+.', '', b[1][0]), b[1][1])])
    return ''.join([x[0][0], '.', x[0][1]]), ''.join([x[1][0], '.', x[1][1]])


def is_number(s):
    """
    Tests if 's' is a number
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def shift(in_list, n=1):
    """
    Shifts a list left (n>0) or right (n<0) in a circular fashion by n positions.
    """
    return in_list[n:] + in_list[:n]


def get_identifier(name, delimiters=('.', '.')):
    """
    Extracts the label (identifier) from an agent name.
    """
    agent_type, identifier = name.split(delimiters[0])[:2]
    if delimiters[1] and delimiters[0] != delimiters[1]:
        identifier = identifier[:-1]
    return agent_type, identifier


def add_identifier(agent_type, id, delimiters=('.', '.')):
    """
    Creates an agent name by joining its type with a label ('id').
    """
    return agent_type + delimiters[0] + id + delimiters[1]


def sort_site_and_bond_lists(mol, s='both'):
    if s == 'both' or s == 'site':
        for st in mol.signature.site_types:
            mol.free_site_list[st].sort(key=lambda x: alphanum_key(x[0]))
            for (i, site) in enumerate(mol.free_site_list[st]):
                mol.free_site_list_idx[st][site] = i  # indices start with 0

    if s == 'both' or s == 'bond':
        for bt in mol.signature.bond_types:
            mol.bond_list[bt].sort(key=lambda x: (alphanum_key(x[0][0]), alphanum_key(x[0][1])))
            for (i, bond) in enumerate(mol.bond_list[bt]):
                mol.bond_list_idx[bt][bond] = i  # indices start with 0


def copy_molecule(X, count=0, id_shift=0, system=None, signature=None, views=None, nav=True, canon=True):
    """
    'Deep copies' a KappaMolecule by using a temporary deep-copy of the agent dictionary to generate
    a new KappaMolecule. This is simpler than attempting to deep-copy the KappaMolecule structure
    directly.
    """
    agents_copy = ujson.loads(ujson.dumps(X.agents))

    if id_shift == 0:
        X_copy = KappaMolecule(agents_copy,
                               count=count,
                               id_shift=id_shift,
                               system=system,
                               sig=signature,
                               s_views=views,
                               nav=nav,
                               canon=canon,
                               init=False)

        # comprehensions are significantly faster than deepcopy()
        X_copy.composition = {k: v for k, v in X.composition.items()}
        X_copy.free_site = {k: v for k, v in X.free_site.items()}
        X_copy.free_site_list = {k: [x for x in X.free_site_list[k]] for k in X.free_site_list}
        X_copy.free_site_list_idx = {k: {x: v for x, v in X.free_site_list_idx[k].items()} for k in X.free_site_list}
        X_copy.bond_type = {k: v for k, v in X.bond_type.items()}
        X_copy.bond_list = {k: [x for x in X.bond_list[k]] for k in X.bond_list}
        X_copy.bond_list_idx = {k: {x: v for x, v in X.bond_list_idx[k].items()} for k in X.bond_list}
        X_copy.bonds = {k: v for k, v in X.bonds.items()}
        X_copy.agent_self_binding = {k: v for k, v in X.agent_self_binding.items()}
        X_copy.unbinding = {k: v for k, v in X.unbinding.items()}
        X_copy.binding = {k: v for k, v in X.binding.items()}
        X_copy.adjacency = {k: [x for x in X.adjacency[k]] for k in X.adjacency}
        X_copy.type_slice = [k for k in X.type_slice]
        X_copy.navigation = {k: v for k, v in X.navigation.items()}
        X_copy.local_views = {k: {x: v for x, v in X.local_views[k].items()} for k in X.local_views}
        X_copy.canonical = X.canonical
        X_copy.sum_formula = X.sum_formula
        X_copy.rarest_type = X.rarest_type
        X_copy.embedding_anchor = X.embedding_anchor
        X_copy.size = len(X_copy.agents)
        # max label
        X_copy.label_counter = int(get_identifier(next(reversed(X_copy.agents)), delimiters=X_copy.id_sep)[1])
    else:
        X_copy = KappaMolecule(agents_copy,
                               count=count,
                               id_shift=id_shift,
                               system=system,
                               sig=signature,
                               s_views=views,
                               nav=nav,
                               canon=canon,
                               init=False)

        X_copy.composition = {k: v for k, v in X.composition.items()}
        X_copy.local_views = {k: {x: v for x, v in X.local_views[k].items()} for k in X.local_views}
        X_copy.canonical = X.canonical
        X_copy.sum_formula = X.sum_formula
        X_copy.rarest_type = X.rarest_type

        X_copy.initialize_light()

    return X_copy


class Kappa:
    """
    A kappa parser based on regular expressions. (Parses correctly correct expressions, but may also parse
    incorrect expressions without crashing... Sorry.)

    Given a molecule or pattern in kappa format, construct its representation as

    agents[name] =
           {
               'iface': { site_name: {'state': state, 'bond': bond label}
               'info': {'id': local id, 'type': agent type, 'sID': SiteSim identifier, 'degree': int n}
               'local_view': local_view
            }

    This is usually passed to KappaMolecule to generate the internal representation used in this package.
    """

    def __init__(self):
        # change these definitions only if you know what you are doing
        self.symbols = r'[_~][a-zA-Z0-9_~+-]+|[a-zA-Z][a-zA-Z0-9_~+-]*'
        self.sID = r'x[0-9]+:'
        self.id_sep = ('.', '.')  # not any of '(?:[_~][a-zA-Z0-9_~+-]+|[a-zA-Z][a-zA-Z0-9_~+-]*)'

        self.label_counter = 0
        self.agents = {}

        # build regex'es

        site_name_re = r'(' + self.symbols + r')'
        internal_state_re = r'({(?:' + self.symbols + r'|[#]' + r')})'
        # binding_id_re = r'(\[(?:.|\d+)\])'
        binding_re = r'(\[(?:.*)\])'  # we still need to parse the bond expression
        # this can be . | # | number | site.agent (a stub)
        self.binding_state_re = \
            re.compile(r'^' + r'(?:\.|_|#|\d+)' + r'|(?:' + self.symbols + r')\.(?:' + self.symbols + r')')
        # using optional lookahead, since the internal state is optional and there is no prescribed order.
        # (gobble up the string with .*)
        self.site_re = \
            re.compile(r'^' + site_name_re + r'(?=.*' + internal_state_re + r')?' + r'(?=.*' + binding_re + r')?.*')
        agent_name_re = r'(' + self.symbols + r')'
        sID_re = r'(' + self.sID + r')?'
        agent_interface_re = r'\(([^()]*)\)'
        # to dissect agents
        self.agent_re = \
            re.compile(r'^' + sID_re + agent_name_re + agent_interface_re + r'$')
        # to find all agents (hence groups are non-capturing), before dissecting them
        self.agents_re = \
            re.compile(r'(?:' + self.sID + r')?(?:' + self.symbols + r')' + r'\([^()]*\)')

    def parser(self, kappa_string, start_label=1):
        self.label_counter = start_label - 1
        self.agents = {}

        expression = re.sub(r'\s+|\t+|\n+', ' ', kappa_string)  # remove line breaks and white matter
        # capture all agents
        match = self.agents_re.findall(expression)

        for agent in match:
            agent_type, identifier, agent_sID, interface = self.parse_agent(agent)
            if agent_sID is None:
                agent_sID = ''
            agent_name = add_identifier(agent_type, identifier, self.id_sep)
            self.agents[agent_name] = {}
            self.agents[agent_name]['iface'] = interface
            self.agents[agent_name]['info'] = {'id': identifier, 'type': agent_type, 'sID': agent_sID, 'degree': -1}
            self.agents[agent_name]['local_view'] = ''

        return self.agents

    def parse_agent(self, agent_expression):
        match = self.agent_re.match(agent_expression)
        if not match:
            sys.exit('Invalid agent declaration <' + agent_expression + '>')
        agent_sID = match.group(1)
        agent_type = match.group(2)
        self.label_counter += 1
        identifier = str(self.label_counter)
        interface = {}

        # parse the agent interface
        iface = match.group(3)
        # since Kappa allows commas or whitespace as separators,
        # swap all commas for spaces and split by whitespace
        sites = iface.replace(',', ' ').split()
        for item in sites:
            try:
                site_name, state, bond = self.parse_site(item)
                interface[site_name] = {'state': state, 'bond': bond}
                # sort interface by site_name
                interface = dict(sorted(interface.items()))
            except ParseFail:
                sys.exit('Could not parse site ' + item + ' in ' + agent_expression)

        return agent_type, identifier, agent_sID, interface

    def parse_site(self, site_expression):
        match = self.site_re.match(site_expression)
        if not match:
            sys.exit('Could not parse site ' + site_expression)
        # return site name, internal state and binding state (without parentheses)
        site_name = match.group(1)
        if match.group(2):  # the modification state; it may be absent, so we need to check
            internal_state = match.group(2)[1:-1]  # remove parens
        else:
            internal_state = '#'  # don't care

        binding_state = '#'  # don't care (absent) by default
        if match.group(3):  # there is an explicit binding state
            binding_expression = match.group(3)[1:-1]  # remove parens
            match = self.binding_state_re.match(binding_expression)  # continue parsing
            if match:
                binding_state = match.group(0)  # either '.' or '#' or number or stub
                # warning: if the site name starts with '_' we have a problem; fix later...
            else:
                sys.exit('Could not parse binding state ' + binding_expression)

        return site_name, internal_state, binding_state

    def decode(self, canon, views):
        """
        Converts a canonical form into a standard kappa expression.
        """

        # alas, reverse the views to index mapping. We use this function rarely, so this is fine;
        # otherwise, reverse elsewhere and store.
        rev_view = {i: k for k, i in views.items()}

        agents = deque()
        # capture all agents
        pos = 0
        for idx in canon.split('.'):
            idx = int(idx)
            if idx > 0:
                # recover local view from index
                agent = rev_view[idx]
                interface = {}
                match2 = self.agent_re.match(agent)
                # agent_sID = match2.group(1)  # there is none in canonical
                agent_type = match2.group(2)
                # the agent interface
                iface = match2.group(3)
                # Since Kappa allows commas or whitespace as separators,
                # swap all commas for spaces and split by whitespace
                sites = iface.replace(',', ' ').split()
                ports = deque()
                for item in sites:
                    site_name, state, bond = self.parse_site(item)
                    interface[site_name] = {'state': state, 'bond': bond}
                    if bond != '.' and bond != '#':
                        other_agent, other_site = bond.split('.')
                        ports.append([agent_type, site_name, other_agent, other_site])
                agents.append({'type': agent_type,
                               'pos': pos,
                               'iface': dict(sorted(interface.items())),
                               'ports': ports,
                               'back': []})
                pos += 1
            else:
                # This is back-edge information.
                # Remember that '-1' in the back edges when we generated the canonical form? Undo it now.
                agents[-1]['back'].append(int(idx) + 1)

        # Do a DFS on the sites (here called ports...)
        stack = deque()
        bond_label = 1
        while agents[0]['ports']:
            stack.append((agents[0]['ports'][-1], agents[0]))
            while stack:
                port1, agent1 = stack.pop()
                agent_type, site_name, other_agent, other_site = port1
                if agent1['back']:
                    # this is the agent at the other end of the back connection
                    agent2 = agents[-agent1['back'][0]]
                    # look for a port in agent1 and the complementary port in agent2
                    for [a, s, a_, s_] in agent1['ports']:
                        if [a_, s_, a, s] in agent2['ports']:
                            # we have a match
                            agent1['iface'][s]['bond'] = bond_label
                            agent2['iface'][s_]['bond'] = bond_label
                            bond_label += 1
                            # annihilate the matching port
                            agent2['ports'].remove([a_, s_, a, s])
                            # also annihilate the matched port
                            agent1['ports'].remove([a, s, a_, s_])
                            agent1['back'] = agent1['back'][1:]
                            if ([a_, s_, a, s], agent2) in stack:
                                stack.remove(([a_, s_, a, s], agent2))
                            if ([a, s, a_, s_], agent1) in stack:
                                stack.remove(([a, s, a_, s_], agent1))
                            if port1 != [a, s, a_, s_]:
                                # put the port back on stack
                                stack.append((port1, agent1))
                            break
                else:
                    for j in range(agent1['pos'] + 1, len(agents)):
                        # find other_agent, other_site, agent_type, site_name
                        agent2 = agents[j]
                        # note the swap (complementarity)
                        if [other_agent, other_site, agent_type, site_name] in agent2['ports']:
                            agent1['iface'][site_name]['bond'] = bond_label
                            agent2['iface'][other_site]['bond'] = bond_label
                            bond_label += 1
                            # annihilate the matching port
                            agent2['ports'].remove([other_agent, other_site, agent_type, site_name])
                            # also annihilate the matched port
                            agent1['ports'].remove(port1)
                            # if not agent2['back']:
                            for q in agent2['ports']:
                                stack.append((q, agent2))
                            break
        ex = ''
        for agent in agents:
            ex += agent['type'] + '('
            for s in agent['iface']:
                ex += s + '[' + str(agent['iface'][s]['bond']) + ']' + '{' + agent['iface'][s]['state'] + '} '
            ex = ex[:-1]
            ex += '), '
        return ex[:-2]


def Canonical2Complex(canonical, views, nav=True, canon=True):
    """
    Wrapper for creating a Kappa molecule from a canonical form.
    """
    _kappa = Kappa()
    expression = _kappa.decode(canonical, views)
    molecule = KappaMolecule(agents=_kappa.parser(expression), s_views={}, nav=nav, canon=canon)
    del _kappa
    return molecule


def KappaComplex(expression, count=0, id_shift=0, system=None, signature=None, views={}, nav=True, canon=True):
    """
    Wrapper for creating a Kappa molecule from an expression.
    """
    # a shortcut for everyday applications
    _kappa = Kappa()
    molecule = KappaMolecule(agents=_kappa.parser(expression),
                             count=count,
                             id_shift=id_shift,
                             system=system,
                             sig=signature,
                             s_views=views,
                             nav=nav,
                             canon=canon)
    del _kappa
    return molecule


class KappaMolecule:
    """
    Constructs the internal representation of a kappa 'molecule'.

    This representation is built from a 'halfway-there' representation (the agent dictionary)
    provided by the Kappa parser or taken from another molecule.

        self.agents[name] =
           {
            'iface': { site_name: {'state': state, 'bond': bond stub}
            'info': {'id': local id, 'type': agent type, 'sID': SiteSim identifier, 'degree': int n}
            'local_view': local_view
            }
        self.adjacency[name] = [ agent1, agent2, ... ]
                self.bonds   = { ( (agent1, site1), (agent2, site2) ) }  # an 'indicator': d[tuple] = 1

        * the interface dictionary of an agent is sorted by site name (needs Python 3.7+)
        * agent names are unique, consisting of type + identifier, eg Axin.42. (including the last dot),
          where the right and left separators (dots by default) are given by self.idsep.
        * self.bonds is a list of unique tuples -- (agent1, site1), (agent2, site2) -- lexicographically
          sorted on agent.
        * bonds are stubs of the form name@site indicating the name of the agent and site
          that anchors the other end of the bond.
          A dictionary has no order by construction, but we can fake an order by iterating through it using
          an ordered list of its keys, whenever order is desired (such as in re-assigning identifiers or
          pretty printing)
        * all types are string, except when otherwise noted.
    """

    def __init__(self,
                 agents=None,
                 count=0,
                 id_shift=0,
                 sig=None,
                 system=None,
                 s_views={},
                 has_views=False,
                 nav=True,
                 canon=True,
                 init=True):

        # change these definitions only if you know what you are doing
        self.bond_sep = '@'
        self.id_sep = ('.', '.')  # not any of '(?:[_~][a-zA-Z0-9_~+-]+|[a-zA-Z][a-zA-Z0-9_~+-]*)'

        # properties of the molecule ---------------------------------------------------
        self.count = count
        self.size = 0
        self.canonical = ''  # the canonicalized expression
        self.composition = {}
        self.sum_formula = ''
        self.rarest_type = ''

        self.system = system
        self.signature = sig

        # we need these to compute reaction propensities
        self.free_site = {}  # number of free sites of a given type
        self.free_site_list = {}  # lists of sites indexed by site type
        self.bond_type = {}  # number of bonds of a given type (w/o labels, as in self.bonds)
        self.bond_list = {}  # lists of bonds indexed by bond type
        self.free_site_list_idx = {}
        self.bond_list_idx = {}
        self.agent_self_binding = {}  # excluded intra-agent binding opportunities; indexed by bond type

        # reaction propensities of molecular species (takes into account the number of instances <count>)
        # these are combinatorial counts multiplied with reaction rate constants, if defined.
        self.unbinding = {}  # propensity of internal bond dissociation; indexed by bond type
        self.binding = {}  # propensity of internal bond formation; indexed by bond type

        # main data structures representing the complex; some redundancy here for convenience

        # we get the 'agents' data structure from the parser
        self.agents = agents

        self.adjacency = {}
        self.bonds = {}
        self.type_slice = []
        self.embedding_anchor = None
        self.navigation = {}
        # flags
        self.nav = nav
        self.canon = canon
        self.has_local_views = has_views
        self.is_pattern = False
        # Local views of the mixture in the context of which expressions are canonicalized
        self.system_views = s_views
        self.local_views = {}

        # auxiliary variables
        self.label_counter = 0  # largest label
        self.next = 0
        self.id_shift = id_shift

        if self.system:  # override so we don't have to set sig and views
            # signature is only used for computing internal reaction propensities
            self.signature = self.system.signature
            self.canon = self.system.canonicalize
            self.nav = False
            if self.system.mixture:
                self.system_views = self.system.mixture.local_views

        # if data are empty, this generates an empty KappaMolecule
        if not self.agents:
            return

        if init:
            self.initialize()

    def initialize(self):
        # size
        self.size = len(self.agents)

        # replace numeric labels of bonds by stubs
        self.stubbify_bonds(id_shift=self.id_shift)

        # max label when labeling is normalized
        self.label_counter = int(get_identifier(next(reversed(self.agents)), delimiters=self.id_sep)[1])
        # get the composition
        self.get_composition()
        self.rarest_type = next(iter(self.composition))

        canonicalize = False
        if self.system:
            if self.size < self.system.size_threshold and self.canon:
                canonicalize = True
            else:
                self.make_adjacency_lists()
                self.canonical = self.system.served_name
                self.system.served_name += 1
        else:
            canonicalize = self.canon

        if canonicalize:
            self.make_adjacency_lists()
            # requires adjacency list
            if not self.has_local_views:
                self.get_local_views()
            self.make_local_view_lists()
            self.canonical = self.canonicalize()

        if self.nav:
            # Get the type lists for site graph matching. This is mostly for offline processing.
            # In simulation, we don't match via graph traversal, but by computing a canonical form.
            for at in self.composition:
                self.type_slice.extend([[name for name in self.agents if self.agents[name]['info']['type'] == at]])
            self.embedding_anchor = self.type_slice[0][0]
            # construct adjacency lists
            if not self.adjacency:
                self.make_adjacency_lists()
            # assemble the navigation list for embeddings
            self.make_navigation_list()

        # calculate reaction propensities
        if self.signature:
            self.internal_reactivity()

    def initialize_light(self):
        # size
        self.size = len(self.agents)

        # replace numeric labels of bonds by stubs
        self.stubbify_bonds(id_shift=self.id_shift)

        # max label when labeling is normalized
        self.label_counter = int(get_identifier(next(reversed(self.agents)), delimiters=self.id_sep)[1])
        # construct adjacency lists
        self.make_adjacency_lists()

        if self.nav:
            # get the type lists for matching
            for at in self.composition:
                self.type_slice.extend([[name for name in self.agents if self.agents[name]['info']['type'] == at]])
            self.embedding_anchor = self.type_slice[0][0]
            # assemble the navigation list for embeddings
            self.make_navigation_list()

        # calculate reaction propensities
        if self.signature:
            self.internal_reactivity()

    def internal_reactivity(self):
        for bt in self.signature.bond_types:
            Xs, Yp = bt
            # internal cycle formation; X.s and Y.p can bind
            if Xs == Yp:  # symmetry correction
                self.binding[bt] = (self.free_site[Xs] * (self.free_site[Xs] - 1)) / 2.
            else:
                self.binding[bt] = self.free_site[Xs] * self.free_site[Yp] - self.agent_self_binding[bt]

            self.unbinding[bt] = self.bond_type[bt]

            if self.system:
                self.binding[bt] *= self.system.rc_bond_formation_intra
                self.unbinding[bt] *= self.system.rc_bond_dissociation[bt]

    def clear_type_lists(self):
        self.free_site_list = {}
        self.free_site_list_idx = {}
        for st in self.signature.site_types:
            self.free_site[st] = 0
            self.free_site_list[st] = []
            self.free_site_list_idx[st] = {}
        self.bond_list = {}
        self.bond_list_idx = {}
        for bt in self.signature.bond_types:
            self.bond_type[bt] = 0
            self.bond_list[bt] = []
            self.bond_list_idx[bt] = {}
            self.agent_self_binding[bt] = 0

    def stubbify_bonds(self, id_shift=0, normalize=True):
        """
        Replaces bond labels with bond stubs.
        """
        if self.signature:
            self.clear_type_lists()

        if not normalize:
            if id_shift == 0:
                self.stubbify_bonds_no_shift()
            else:
                self.stubbify_bonds_with_shift(id_shift=id_shift)
        else:
            self.stubbify_bonds_with_shift(remap=self.normalize_ids(id_shift=id_shift))

        # sort_site_and_bond_lists(self)

    def stubbify_bonds_no_shift(self):
        """
        Replaces numeric bond labels with unique bond stubs.
        For example, A.14.(b[2]), Z.3.(j[2]) becomes A.14.(b[Z.3.@j]), Z.3.(j[A.14.@b]).

        generates:
            self.bonds
            self.bond_type
            self.bond_type_list
            self.free_site
            self.free_site_list
            self.agent_self_binding
        """
        # Note: If we are dealing with an object that contains a bond pattern, the degree of a node has no meaning.
        self.bonds = {}
        bonds = {}
        stubs = []
        for name in self.agents:
            degree = 0
            agent_free_site_types = set()
            for site in self.agents[name]['iface']:
                link = self.agents[name]['iface'][site]['bond']
                if link != '.':
                    if is_number(link):
                        degree += 1
                        if link in bonds:
                            [(name1, site1)] = bonds[link]
                            # stubbify
                            self.agents[name1]['iface'][site1]['bond'] = ''.join([name, self.bond_sep, site])
                            self.agents[name]['iface'][site]['bond'] = ''.join([name1, self.bond_sep, site1])
                        else:
                            bonds[link] = [(name, site)]
                            continue
                    # this occurs when we created the molecule with an already 'stubbified' agent dictionary
                    elif self.bond_sep in link:
                        degree += 1
                        name1, site1 = link.split(self.bond_sep)
                        complement = self.agents[name1]['iface'][site1]['bond']
                        if complement not in stubs:
                            stubs += [link]
                            continue
                    else:
                        # bond state is a ghost, or '_', or '#'
                        # degree = -1  # reset and flag, just in case
                        self.is_pattern = True
                        continue
                    # This purpose of this section is to
                    #  (1) identify the bonds and store them as keys in self.bonds
                    #  (2) count the bond types (used to calculate reaction propensities).
                    # Note: we can get here only from the case in which we have already seen the partner
                    # of the bond of agent 'name' at site 'site'. That partner is 'name1' at site 'site1'.
                    # We standardize the bond by sorting.
                    b = sorted([(name1, site1), (name, site)], key=lambda x: (alphanum_key(x[0]), alphanum_key(x[1])))
                    b = tuple(b)
                    # collect unique bonds
                    self.bonds[b] = 1  # just an indicator; we are collecting unique keys (bonds)
                    # count the bond *types* (to compute reactivity); here labels don't matter
                    if self.signature:
                        bt = bond2type(b)
                        self.bond_list[bt].append(b)
                        self.bond_list_idx[bt][b] = self.bond_type[bt]  # (ab)used as a counter
                        self.bond_type[bt] += 1
                else:
                    # The site is not bound.
                    # Accumulate the free sites in the molecule.
                    if self.signature:
                        # count the free sites (to compute reaction propensities)
                        st = ''.join([self.agents[name]['info']['type'], '.', site])
                        agent_free_site_types.add(st)  # a local agent-specific set, used below
                        self.free_site_list[st].append((name, site))
                        self.free_site_list_idx[st][(name, site)] = self.free_site[st]  # (ab)used as a counter
                        self.free_site[st] += 1
            if self.signature:
                # Accumulate the intra-agent binding opportunities within the whole molecule.
                # This will be used to correct the intra-molecular bond formation propensity
                # (by removing agent self-binding).
                for bt in self.signature.bond_types:
                    st1, st2 = bt
                    # The case st1 == st2 cannot arise within a single Kappa agent
                    if st1 != st2:
                        if st1 in agent_free_site_types and st2 in agent_free_site_types:
                            self.agent_self_binding[bt] += 1

            self.agents[name]['info']['degree'] = degree

    def stubbify_bonds_with_shift(self, id_shift=0, remap=None):
        """
        Replaces numeric bond labels with unique bond stubs much like stubbify_bonds_no_shift(),
        but also executes a label remapping and label shift. Absent a remap and an id_shift,
        it amounts to stubbify_bonds_no_shift(). However, since it creates a new agent dictionary,
        it is preferable to use stubbify_bonds_no_shift() in that case, especially if
        complexes are very large. (We could merge the two functions at the cost of a few conditionals.)

        The label shift is needed when we connect molecules into a larger molecules.
        The fusion requires shifting the agent identifiers of one of the molecules by the number
        of agents contained in the other.
        """
        self.bonds = {}
        bonds = {}
        stubs = []
        new_agents = {}
        remapping = remap
        if not remap:
            # identity
            remapping = [str(i + id_shift) for i in range(1, len(self.agents) + 1)]

        for name in self.agents:
            #
            type1 = self.agents[name]['info']['type']
            sID = self.agents[name]['info']['sID']
            new_id = remapping[self.agents[name]['info']['id']]
            # new_id = str(int(self.agents[name]['info']['id']) + id_shift)
            new_name = add_identifier(self.agents[name]['info']['type'], new_id)
            new_agents[new_name] = {}
            new_agents[new_name]['iface'] = {}
            new_agents[new_name]['info'] = {'id': new_id, 'type': type1, 'sID': sID, 'degree': 0}
            new_agents[new_name]['local_view'] = self.agents[name]['local_view']  # shift-independent
            degree = 0
            agent_free_site_types = []
            for site in self.agents[name]['iface']:
                new_agents[new_name]['iface'][site] = {}
                new_agents[new_name]['iface'][site]['state'] = self.agents[name]['iface'][site]['state']
                # may be overwritten below
                new_agents[new_name]['iface'][site]['bond'] = self.agents[name]['iface'][site]['bond']
                link = self.agents[name]['iface'][site]['bond']
                if link != '.':
                    if is_number(link):
                        degree += 1
                        if link in bonds:
                            (name1, site1) = bonds[link]
                            # stubbify
                            new_id1 = remapping[self.agents[name1]['info']['id']]
                            # new_id1 = str(int(self.agents[name1]['info']['id']) + id_shift)
                            new_name1 = add_identifier(self.agents[name1]['info']['type'], new_id1)
                            new_agents[new_name1]['iface'][site1]['bond'] = ''.join([new_name, self.bond_sep, site])
                            new_agents[new_name]['iface'][site]['bond'] = ''.join([new_name1, self.bond_sep, site1])
                        else:
                            bonds[link] = (name, site)
                            continue
                    # this occurs when we created the molecule with an already 'stubbified' agent dictionary
                    elif self.bond_sep in link:
                        degree += 1
                        name1, site1 = link.split(self.bond_sep)
                        complement = self.agents[name1]['iface'][site1]['bond']
                        if complement in stubs:
                            new_id1 = remapping[self.agents[name1]['info']['id']]
                            # new_id1 = str(int(self.agents[name1]['info']['id']) + id_shift)
                            new_name1 = add_identifier(self.agents[name1]['info']['type'], new_id1)
                            new_agents[new_name1]['iface'][site1]['bond'] = ''.join([new_name, self.bond_sep, site])
                            new_agents[new_name]['iface'][site]['bond'] = ''.join([new_name1, self.bond_sep, site1])
                        else:
                            stubs += [link]
                            continue
                    else:
                        # bond state is a ghost, or '_', or '#'
                        # degree = -1  # reset and flag, just in case
                        self.is_pattern = True
                        continue
                    n1 = self.agents[name1]['info']['type']
                    n = self.agents[name]['info']['type']
                    (t1, l1, s1), (t2, l2, s2) = sorted([(n1, int(new_id1), site1), (n, int(new_id), site)])
                    b = (add_identifier(t1, str(l1)), s1), (add_identifier(t2, str(l2)), s2)
                    # collect unique bonds
                    self.bonds[b] = 1  # just an indicator; we are collecting unique keys (bonds)
                    # count the bond *types* (to compute reactivity); here labels don't matter
                    if self.signature:
                        (t1, s1), (t2, s2) = sorted([(n1, site1), (n, site)])
                        bt = (''.join([t1, '.', s1]), ''.join([t2, '.', s2]))
                        self.bond_list[bt].append(b)
                        self.bond_list_idx[bt][b] = self.bond_type[bt]  # (ab)used as a counter
                        self.bond_type[bt] += 1
                else:
                    if self.signature:
                        # count the free sites for the whole molecule (to compute reactivity)
                        st = ''.join([self.agents[name]['info']['type'], '.', site])
                        agent_free_site_types += [st]
                        self.free_site_list[st].append((new_name, site))
                        self.free_site_list_idx[st][(new_name, site)] = self.free_site[st]  # (ab)used as a counter
                        self.free_site[st] += 1
            if self.signature:
                for bt in self.signature.bond_types:
                    st1, st2 = bt
                    if st1 != st2:
                        if st1 in agent_free_site_types and st2 in agent_free_site_types:
                            self.agent_self_binding[bt] += 1

            new_agents[new_name]['info']['degree'] = degree

        self.agents = new_agents

    def get_local_views(self):
        """
        Obtain the lexically ordered local view at each agent.
        """
        # get the last used system-wide index of local views encountered thus far
        if not self.system_views:
            running_id = 0
        else:
            running_id = self.system_views[next(reversed(self.system_views))]

        # get the "local views"
        self.local_views = {}
        for name in self.agents:
            lv = []
            iface = self.agents[name]['iface']
            for s in iface:
                view = ''
                b = iface[s]['bond']
                if b != '.' and b != '#':
                    other_name, other_s = b.split(self.bond_sep)
                    other_type = self.agents[other_name]['info']['type']
                    view += f'[{other_type}.{other_s}]'
                else:
                    view += f'[{b}]'
                # skip the state in this specific context
                # view += '{' + f"{iface[s]['state']}" + '}'
                lv.append((s, view))
            local_view = ''
            for site_view in [f"{s}{view} " for (s, view) in sorted(lv)]:
                local_view += site_view
            # this is the local view of agent 'name'
            l_view = self.agents[name]['info']['type'] + '(' + local_view[:-1] + ')'

            self.agents[name]['local_view'] = l_view

            # update the system views
            if l_view not in self.system_views:
                running_id += 1
                self.system_views[l_view] = running_id

    def make_local_view_lists(self):
        self.local_views = {}
        for name in self.agents:
            # make lists of agents with the same local view
            lv = self.agents[name]['local_view']
            if lv in self.local_views:
                self.local_views[lv][name] = 1  # just an indicator dict for fast search and deletion
            else:
                self.local_views[lv] = {name: 1}

    def canonicalize(self):
        """
        Canonicalize the kappa expression.
        """
        if not self.system_views:
            return ''
        # get the local view with the smallest index in the _system_ (!)
        _, mlv = min([(self.system_views[lv], lv) for lv in self.local_views])
        # This is the list of local nodes with that view
        node_list = list(self.local_views[mlv].keys())

        canonic = []
        for node in node_list:
            canonic.append(self.traverse(node))
        # return the lexicographically smallest
        return min(canonic)

    def traverse(self, node):
        """
        Makes a DFS traversal and identifies back edges, then constructs the canonical form as
        a sequence of integers that are indices to the local views. Negative integers are back edges
        to the position in the sequence (after adding 1 [used to avoid zero]). Assumes an undirected graph.
        """
        discovered = set()
        spanning = set()
        cycle_edges = defaultdict(int)
        traversal = []
        traversal_index = {}
        idx = 0
        parent = {node: node}
        stack = deque()
        stack.append((node, node))
        while stack:
            # pop the 'current' node and the node 'p' that preceded it when 'current' was placed on the stack
            current, p = stack.pop()
            if current not in discovered:
                spanning.add(tuple(sorted((parent[current], current))))  # sort, since graph is undirected
                # the list associated with 'current' is to later hold the back edges
                traversal.append((current, []))  # the sequence is important, obviously
                # we use a dictionary to store the positions in 'traversal' for fast retrieval below
                traversal_index[current] = idx
                idx += 1
                discovered.add(current)  # for faster search (or else we could use 'traversal')
                for neighbor in self.adjacency[current]:
                    if neighbor not in discovered:
                        parent[neighbor] = current
                        stack.append((neighbor, current))
            else:
                e = tuple(sorted((p, current)))
                # 'cycle_edges' is a counter. This handles multi-graphs, meaning graphs in which two nodes
                # can have multiple edges between them.
                cycle_edges[e] += 1
        # --------------------------------------------------------------------------------
        # insert the back edges
        for e in cycle_edges:
            for _ in range(cycle_edges[e]):
                if traversal_index[e[0]] > traversal_index[e[1]]:
                    # back edge to e[1]
                    (current, back_list) = traversal[traversal_index[e[0]]]
                    traversal[traversal_index[e[0]]] = (current, back_list + [e[1]])
                else:
                    # back edge to e[0]
                    (current, back_list) = traversal[traversal_index[e[1]]]
                    traversal[traversal_index[e[1]]] = (current, back_list + [e[0]])

        # convert into a proper list and re-index
        canonic = []
        idx = 0
        for (n, back_list) in traversal:
            canonic.append(str(self.system_views[self.agents[n]['local_view']]))
            traversal_index[n] = idx
            idx += 1
            if back_list:
                for i in back_list:  # sorting is necessary to achieve a deterministic result
                    # Make sure you get a negative number, i.e. not 0, as 0 might conflict with a
                    # local-view index of 0... Although our indices start with 1 (and thus no conflict), we
                    # subtract 1 to get a visual minus sign. Don't forget to +1 to get back the index...
                    canonic.append(str(-traversal_index[i] - 1))

        return '.'.join(canonic)

    def make_adjacency_lists(self):
        """
        Construct adjacency lists for each agent
        """
        self.adjacency = {}
        for name1 in self.agents:
            iface = self.agents[name1]['iface']
            adjacency = [iface[s1]['bond'].split(self.bond_sep)[0] for s1 in iface if iface[s1]['bond'] != '.']
            self.adjacency[name1] = adjacency

    def make_navigation_list(self):
        # self.navigation[(a1, a2)] contains a site of a1 that anchors a bond to a2
        # (For the purpose of this array, we don't care about multiple bonds between the same agents.)
        # This is similar to self.bonds, but organized as a dictionary for convenience.
        self.navigation = {}
        for (a1, s1), (a2, s2) in self.bonds:  # names a1 and a2 in bonds have 'id' attached
            self.navigation[(a1, a2)] = s1
            self.navigation[(a2, a1)] = s2

    def shift_ids(self, id_shift=0):
        remapping = {}
        for name in self.agents:
            ident = self.agents[name]['info']['id']
            remapping[ident] = str(int(ident) + id_shift)
        return remapping

    def normalize_ids(self, id_shift=0):
        remapping = {}
        i = 1
        for name in self.agents:
            remapping[self.agents[name]['info']['id']] = str(i + id_shift)
            i += 1
        return remapping

    def randomize_ids(self):
        l = [i for i in range(1, len(self.agents) + 1)]
        random.shuffle(l)
        # random.Random(42).shuffle(l)
        remapping = {}
        i = 0
        for name in self.agents:
            remapping[self.agents[name]['info']['id']] = str(l[i])
            i += 1
        return remapping

    def remap_ids(self, remapping):
        """
        (Re)assigns agent labels (identifiers) using the map 'remapping'.
        """
        self.bonds = {}  # reset

        # apply permutation
        renamed = {}
        for name1 in self.agents:
            id1 = self.agents[name1]['info']['id']
            type1 = self.agents[name1]['info']['type']
            sID = self.agents[name1]['info']['sID']
            new_id1 = remapping[id1]
            new_name1 = add_identifier(type1, new_id1)
            renamed[new_name1] = {}
            renamed[new_name1]['iface'] = {}
            renamed[new_name1]['info'] = {'id': new_id1, 'type': type1, 'sID': sID,
                                          'degree': self.agents[name1]['info']['degree']}
            renamed[new_name1]['local_view'] = ''
            for site1 in self.agents[name1]['iface']:
                renamed[new_name1]['iface'][site1] = {}
                renamed[new_name1]['iface'][site1]['state'] = self.agents[name1]['iface'][site1]['state']
                if self.bond_sep in self.agents[name1]['iface'][site1]['bond']:
                    # agent2 is name.old_id
                    agent2, site2 = self.agents[name1]['iface'][site1]['bond'].split(self.bond_sep)
                    type2, id2 = get_identifier(agent2)
                    new_id2 = remapping[id2]
                    new_name2 = add_identifier(type2, new_id2)
                    renamed[new_name1]['iface'][site1]['bond'] = new_name2 + self.bond_sep + site2
                    # sort
                    (t1, l1, s1), (t2, l2, s2) = sorted([(type1, int(new_id1), site1), (type2, int(new_id2), site2)])
                    b = (add_identifier(t1, str(l1)), s1), (add_identifier(t2, str(l2)), s2)
                    # collect unique bonds
                    if b not in self.bonds:
                        self.bonds[b] = 1  # just an indicator; we are collecting unique keys (bonds)
                else:
                    renamed[new_name1]['iface'][site1]['bond'] = self.agents[name1]['iface'][site1]['bond']

        self.agents = renamed

        if self.signature:
            new_bond_list = {}
            new_bond_list_idx = {}
            for bt in self.signature.bond_types:
                new_list = []
                new_list_idx = {}
                for (b1, b2) in self.bond_list[bt]:
                    n1, s1 = b1
                    n2, s2 = b2
                    t1, id1 = get_identifier(n1)
                    t2, id2 = get_identifier(n2)
                    new_id1 = remapping[id1]
                    new_id2 = remapping[id2]
                    (t1, l1, s1), (t2, l2, s2) = sorted([(t1, int(new_id1), s1), (t2, int(new_id2), s2)])
                    b = (add_identifier(t1, str(l1)), s1), (add_identifier(t2, str(l2)), s2)
                    new_list.append(b)
                    new_list_idx[b] = len(new_list) - 1
                new_bond_list[bt] = new_list
                new_bond_list_idx[bt] = new_list_idx
            self.bond_list = new_bond_list
            self.bond_list_idx = new_bond_list_idx

            new_free_site_list = {}
            new_free_site_list_idx = {}
            for st in self.signature.site_types:
                new_list = []
                new_list_idx = {}
                for (name, site) in self.free_site_list[st]:
                    t1, id1 = get_identifier(name)
                    new_id1 = remapping[id1]
                    new_name = add_identifier(t1, str(new_id1))
                    new_list.append((new_name, site))
                    new_list_idx[(new_name, site)] = len(new_list) - 1
                new_free_site_list[st] = new_list
                new_free_site_list_idx[st] = new_list_idx
            self.free_site_list = new_free_site_list
            self.free_site_list_idx = new_free_site_list_idx

    def remap(self, change='none', id_shift=0):
        """
        A wrapper for remap_ids() -- of use for external calls; not used in setting up the object.
        'change' = {'none', 'normalize', 'randomize'} directs the construction of the remapping map.
        Identifiers are shifted by 'id_shift'.
        """
        if change == 'normalize':
            self.remap_ids(self.normalize_ids())
        elif change == 'randomize':
            self.remap_ids(self.randomize_ids())
        else:  # check if shift
            if id_shift > 0:
                self.remap_ids(self.shift_ids(id_shift=id_shift))

        # construct adjacency lists
        self.make_adjacency_lists()

        if self.nav:
            # get the type lists for matching
            self.type_slice = []
            for at in self.composition:
                self.type_slice.extend([[name for name in self.agents if self.agents[name]['info']['type'] == at]])
            self.embedding_anchor = self.type_slice[0][0]
            # assemble the navigation list for embeddings
            self.make_navigation_list()

        if self.canon:
            self.get_local_views()
            self.make_local_view_lists()
            self.canonical = self.canonicalize()

    def get_composition(self):
        """
        Get the 'sum formula' of a complex. Agents are ordered in increasing abundance within the complex.
        """
        comp = {}
        for a in self.agents:
            a_type = self.agents[a]['info']['type']
            if a_type in comp:
                comp[a_type] += 1
            else:
                comp[a_type] = 1

        # sort the dict by value and then key:
        self.composition = {k: v for k, v in sorted(comp.items(), key=lambda item: (item[1], item[0]))}

        self.sum_formula = ''
        for a_type in self.composition:
            self.sum_formula += (a_type + '{' + str(self.composition[a_type]) + '}')

    def is_multigraph(self):
        """
        Test if the set of bonds implies a multi-graph.
        """
        s = set()
        for (a1, s1), (a2, s2) in self.bonds:
            if (a1, a2) in s:
                return True
            else:
                s.add((a1, a2))
        return False

    def nodes(self):
        """
        Emulates the networkx G.nodes() method returning a list of node names.
        """
        return [k for k in self.agents]

    def order(self):
        """
        Works like __len__. For compatibility with networkx representation.
        """
        return self.size

    def degree(self):
        """
        Emulates networkx G.degree(), returning a list of (node, degree) pairs
        """
        l = []
        for name in self.agents:
            if self.agents[name]['info']['degree'] == -1:
                return []
            l += [(name, self.agents[name]['info']['degree'])]
        return l

    def kappa_expression(self, label=False):
        """
        Converts the internal representation of a kapa molecule into a kappa string
        """
        # If we are dealing with an object that contains a bond pattern, the degree of a node has no meaning.
        # The degree is used only for VF2 isomorphism checking, but not for pattern embeddings.
        i = 1
        num = {}
        s = ''
        for name in self.agents:
            s += self.agents[name]["info"]["sID"]
            if label:
                s += f'{name}('
            else:
                s += f'{self.agents[name]["info"]["type"]}('
            for site in self.agents[name]['iface']:
                s += f'{site}'
                state = self.agents[name]['iface'][site]['state']
                if state != '#':
                    s += '{' + f'{state}' + '}'
                link = self.agents[name]['iface'][site]['bond']
                if link == '.':
                    s += '[.] '
                elif link == '#':
                    s += '[#] '
                    # s += ''
                elif self.bond_sep in link:
                    ag, ste = link.split(self.bond_sep)
                    if (name, site) in num:
                        s += f'[{num[(name, site)]}] '
                    else:
                        num[(ag, ste)] = i
                        s += f'[{i}] '
                        i += 1
                else:
                    s += f'[{link}] '
            if not self.agents[name]['iface']:
                s = s + '), '
            else:
                s = s[:-1] + '), '
        return s[:-2]

    def summary(self, internal=False, show_bonds=False, reactivity=False, db_level=0, pp_width=40):
        """
        Prints summary of the molecule at various levels of detail.
        """
        info = '\n'
        info += f"{''.ljust(70, '-')}\n"
        n_free_sites = 0
        if not self.signature:
            info += 'Warning: no signature. Counts presume a signature.\n'
        for st in self.free_site:
            n_free_sites += self.free_site[st]
        info += f'[count: {self.count}] {self.size} agents, {len(self.bonds)} bonds, and {n_free_sites} free sites\n'
        info += f'composition: {self.sum_formula}\n'
        if self.is_pattern:
            info += f'expression is a pattern !\n'
        info += self.show(internal=internal, wrap=200) + '\n'
        info += f"{''.ljust(70, '-')}\n"
        if show_bonds:
            info += self.report_bond_types_and_free_sites(db_level=db_level, pp_width=3)
        if reactivity:
            info += self.report_internal_reaction_propensities(pp_width=pp_width)
        return info

    def report_bond_types_and_free_sites(self, db_level=0, pp_width=40):
        """
        Prints the bond types and free site types of the molecule.
        """
        info = ''
        if self.signature:
            s = f'{len(self.bond_type)} bond types:'
            info += f'{s:>{pp_width}}\n'
            for bt in self.bond_type:
                if self.bond_type[bt] != 0:
                    s1, s2 = bt
                    b = f"{s1}-{''.join([s2.split('.')[1], '.', s2.split('.')[0]])}"
                    info += f'{"":>{pp_width}} {b}: {self.bond_type[bt]}\n'
                    temp = f'{"":>{pp_width}} {b}: {self.bond_list[bt]}'
                    temp = pprint.pformat(temp, indent=0, width=200, compact=False)
                    info += temp[1:-1].replace('"', '') + '\n'
            s = f'{len(self.free_site)} free site types:'
            info += f'{s:>{pp_width}}\n'
            for st in self.free_site:
                if self.free_site[st] != 0:
                    info += f'{"":>{pp_width}} {st}: {self.free_site[st]}\n'
                    info += f'{"":>{pp_width}} {st}: {self.free_site_list[st]}\n'
        else:
            info = 'Warning: no signature. (Use db_level 2 for a full list of bonds.)\n'

        if db_level == 2:
            info += f"\n"
            s = f'list of bonds (random order):'
            info += f'{s:>{pp_width}}\n'
            for (a1, s1), (a2, s2) in self.bonds:
                b1 = ''.join([a1, s1])
                b2 = ''.join([a2, s2])
                info += f'{"":>{pp_width}} {b1}<->{b2}\n'
        return info

    def report_internal_reaction_propensities(self, pp_width=40):
        """
        Prints the intra-molecular binding propensities of the molecule.
        """
        if self.signature:
            form = '1.5E'
            info = "\n"
            s = f'reaction propensities'
            info += f'{s:>{pp_width}}\n'
            info += f'{"dissociation per instance":>{pp_width}}:\n'
            for bt in self.unbinding:
                s1, s2 = bt
                b = f"{s1}-{''.join([s2.split('.')[1], '.', s2.split('.')[0]])}"
                info += f'{b:>{pp_width}}: {self.unbinding[bt]:{form}}\n'
            info += f'{"binding per instance":>{pp_width}}\n'
            for bt in self.binding:
                s1, s2 = bt
                b = f"{s1}-{''.join([s2.split('.')[1], '.', s2.split('.')[0]])}"
                info += f'{b:>{pp_width}}: {self.binding[bt]:{form}}\n'
            info += f'{"excluded agent self-binding counts per instance":>{pp_width}}\n'
            for bt in self.agent_self_binding:
                if self.agent_self_binding[bt] != 0:
                    s1, s2 = bt
                    b = f"{s1}-{''.join([s2.split('.')[1], '.', s2.split('.')[0]])}"
                    info += f'{b:>{pp_width}}: {self.agent_self_binding[bt]}\n'
        else:
            info = f'\n{"Warning: no signature":>30}'

        return info

    def show(self, internal=False, label=False, wrap=-1):
        """
        Prints the internal representation.
        If internal=False, print the standard kappa expression.
        """
        info = ''
        if internal:
            for name in self.agents:
                interface = ''
                iface = self.agents[name]['iface']
                for s in iface:
                    interface += s + '{' + iface[s]['state'] + '}' + '[' + iface[s]['bond'] + '] '
                if not self.is_pattern:
                    info += f"[d: {self.agents[name]['info']['degree']}] "
                info += self.agents[name]['info']['sID'] + name + '(' + interface[:-1] + ')\n'
            return info[:-1]
        else:
            if wrap > 0:
                info = pprint.pformat(self.kappa_expression(label=label), indent=0, width=wrap, compact=False)
                return info[1:-1].replace("'", "")
            else:
                return self.kappa_expression(label=label)

    def __str__(self):
        return self.kappa_expression()

    def __iter__(self):
        return iter(self.agents)

    def __len__(self):
        return self.size

    def __getitem__(self, name):
        """
        Makes C[name] return the list of neighbors of node name; emulates the adjacency view of networkx
        """
        return self.adjacency[name]

    def show_local_views(self):
        local_view = {}
        info = ''
        for n in self.agents:
            if self.agents[n]['local_view'] in local_view:
                local_view[self.agents[n]['local_view']] += 1
            else:
                local_view[self.agents[n]['local_view']] = 1
        for lv in local_view:
            info += f"{local_view[lv]:3d}  {lv}\n"
        return info


# -------------------------------------------------------------------------------------------


if __name__ == '__main__':
    import kasnap as snap
    import kamatch

    # parser
    kappa = Kappa()
    # # a simple Kappa string
    # s1 = ' A(o[1], p[2] t{p}[3]), B(x[1] y[2] z[.]), C(w[3], y[z.B])'
    # agents = kappa.parser(s1)
    # c = KappaMolecule(agents, count=175)
    # out = c.kappa_expression()
    # print(f"expression:\n{out}")
    # out = c.show(internal=True)
    # print(f"internal representation:\n{out}")

    print("--------------")

    s2 = " x222667:P(a1[.] a2[.] a3[.] d[1]), x251065:P(a1[.] a2[.] a3[.] d[1])"
    agents = kappa.parser(s2)
    c = KappaMolecule(agents)
    out = c.kappa_expression()
    print(f"expression:\n{out}")
    out = c.show(internal=True)
    print(f"internal representation:\n{out}")

    print("--------------")

    print("Reading snapshot")
    snapshot = snap.SnapShot(file='TestData/snap__1784.ka')
    print("Done reading")
    print(snapshot.snap_report())

    # print("--------------")
    #
    SGM = kamatch.SiteGraphMatcher()
    #
    # for molecule in snapshot.complexes:
    #     canonical = molecule.canonical
    #     print(canonical)
    #     molecule_ = KappaComplex(kappa.decode(canonical, snapshot.local_views))
    #     print(f'decoded is isomorphic to original: {SGM.isomorphic(molecule, molecule_)}')

    print("--------------")

    print("Reading snapshot")
    snapshot = snap.SnapShot(file='TestData/snap__1773.ka')
    print("Done reading")
    print(snapshot.snap_report())

    example = snapshot.complexes[0]
    print(example.show(label=True, wrap=100))
    print("same thing")
    s = KappaMolecule(example.agents)
    print(s.show(label=True, wrap=100))

    print("> remapped______________")
    example.remap(change='randomize', id_shift=100)
    print(example.show(label=True, wrap=100))

    s = KappaMolecule(example.agents)
    print("agent ids have been normalized")
    print(s.show(label=True, wrap=100))

    print(s.summary(internal=True, show_bonds=True, reactivity=True))

    # import kasig as sig

    data = 'A(a1[1] a2[2] a3[3] c[8]), A(a1[1] a2[2] a3[3] c[4]) A(a1[5] a2[6] a3[7] c[4]), A(a1[5] a2[6] a3[7] c[8])'
    x1 = KappaComplex(data)
    print(x1.kappa_expression())
    print(x1.canonical)
    out = kappa.decode(x1.canonical, x1.system_views)
    print(out)
    print(f'decoded is isomorphic to original: {SGM.isomorphic(KappaComplex(out), x1)}')
