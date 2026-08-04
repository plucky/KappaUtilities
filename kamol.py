# Walter Fontana at 04/28/23
# Reorganized by Claude Caude (Sonnet5 High) 08/04/2026
"""
This module defines the internal representation of a kappa molecule. It implements the parser, and a
 number of utilities, such as copying molecules
"""
import re
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


def add_identifier(agent_type, identifier, delimiters=('.', '.')):
    """
    Creates an agent name by joining its type with a label ('identifier').
    """
    return agent_type + delimiters[0] + identifier + delimiters[1]


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


def copy_molecule(X, count=0, id_shift=0, system=None, signature=None, local_view_index=None, nav=True, canon=True):
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
                               signature=signature,
                               local_view_index=local_view_index,
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
                               signature=signature,
                               local_view_index=local_view_index,
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

    def parse(self, kappa_string, start_label=1):
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
            raise ParseFail(f'Invalid agent declaration <{agent_expression}>')
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
            except ParseFail:
                raise ParseFail(f'Could not parse site {item} in {agent_expression}')
            interface[site_name] = {'state': state, 'bond': bond}
            # sort interface by site_name
            interface = dict(sorted(interface.items()))

        return agent_type, identifier, agent_sID, interface

    def parse_site(self, site_expression):
        match = self.site_re.match(site_expression)
        if not match:
            raise ParseFail(f'Could not parse site {site_expression}')
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
                raise ParseFail(f'Could not parse binding state {binding_expression}')

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


def canonical_to_expression(canonical, local_view_index, nav=True, canon=True):
    """
    Wrapper for creating a KappaExpression from a canonical form.
    """
    parser = Kappa()
    k_expression = parser.decode(canonical, local_view_index)
    return KappaExpression(agents=parser.parse(k_expression), nav=nav, canon=canon)


def kappa_to_expression(k_expression, id_shift=0, nav=True, canon=False):
    """
    Wrapper for creating a KappaExpression from a kappa string. A shortcut for everyday applications.
    """
    parser = Kappa()
    return KappaExpression(agents=parser.parse(k_expression), id_shift=id_shift, nav=nav, canon=canon)


def kappa_to_molecule(k_expression, count=0, id_shift=0, system=None, signature=None, local_view_index=None,
                      nav=True, canon=True):
    """
    Wrapper for creating a KappaMolecule from a kappa string. A shortcut for everyday applications.
    """
    parser = Kappa()
    return KappaMolecule(agents=parser.parse(k_expression), count=count, id_shift=id_shift,
                          system=system, signature=signature, local_view_index=local_view_index,
                          nav=nav, canon=canon)


class KappaExpression:
    """
    Constructs the internal representation of a kappa expression.

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

    def __init__(self, agents=None, id_shift=0, nav=False, canon=False, init=True):

        # change these definitions only if you know what you are doing
        self.bond_sep = '@'
        self.id_sep = ('.', '.')  # not any of '(?:[_~][a-zA-Z0-9_~+-]+|[a-zA-Z][a-zA-Z0-9_~+-]*)'

        # properties of the complex ---------------------------------------------------
        self.size = 0
        self.n_free_sites = 0
        self.canonical = ''  # the canonicalized expression
        self.composition = {}
        self.sum_formula = ''
        self.rarest_type = ''

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
        self.is_pattern = False
        # Local views for canonicalization
        self.local_views = {}
        self.local_view_index = {}

        # auxiliary variables
        self.label_counter = 0  # largest label
        self.next = 0
        self.id_shift = id_shift

        if init:
            self.initialize_expression(canon=self.canon, nav=self.nav)

    def initialize_expression(self, canon=False, nav=False):
        # size
        self.size = len(self.agents)

        # replace numeric labels of bonds by stubs
        self.stubbify_bonds(id_shift=self.id_shift)

        # max label when labeling is normalized
        self.label_counter = int(get_identifier(next(reversed(self.agents)), delimiters=self.id_sep)[1])
        # get the composition
        self.get_composition()
        self.rarest_type = next(iter(self.composition))

        if canon and not self.is_pattern:
            self.make_adjacency_lists()
            # requires adjacency list
            self.get_local_views()
            self.make_local_view_lists()
            self.canonical = self.canonicalize()

        if nav:
            # Get the type lists for site graph matching. This is mostly for offline processing.
            # In simulation, we don't match via graph traversal, but by computing a canonical form.
            self._setup_navigation()

    def _setup_navigation(self):
        """
        Builds per-type agent-name lists (self.type_slice) and the site-graph navigation table
        (self.navigation) used for pattern-matching embeddings.
        """
        self.type_slice = [[name for name in self.agents if self.agents[name]['info']['type'] == at]
                           for at in self.composition]
        self.embedding_anchor = self.type_slice[0][0]
        if not self.adjacency:
            self.make_adjacency_lists()
        self.make_navigation_list()

    def stubbify_bonds(self, id_shift=0, normalize=True):
        """
        Replaces numeric bond labels in agent interfaces with unique bond stubs of the form
        'name@site', which is the bond representation used throughout this module.
        """
        self.n_free_sites = 0
        self.is_pattern = False

        if not normalize:
            if id_shift == 0:
                self._stubbify_inplace()
            else:
                self._stubbify_rebuild(id_shift=id_shift)
        else:
            self._stubbify_rebuild(remapping=self.normalize_ids(id_shift=id_shift))

    def _stubbify_inplace(self):
        """
        Replaces numeric bond labels with unique bond stubs, in place.
        For example, A.14.(b[2]), Z.3.(j[2]) becomes A.14.(b[Z.3.@j]), Z.3.(j[A.14.@b]).
        """
        # Note: If we are dealing with an object that contains a bond pattern, the degree of a node has no meaning.
        self.bonds = {}
        pending = {}  # numeric bond label -> (name, site) awaiting its partner
        confirmed_stubs = []  # bond stubs seen once, awaiting confirmation from their partner site
        for name in self.agents:
            degree = 0
            for site in self.agents[name]['iface']:
                link = self.agents[name]['iface'][site]['bond']
                if link == '.':
                    self.n_free_sites += 1
                    continue
                if is_number(link):
                    degree += 1
                    if link in pending:
                        [(name1, site1)] = pending[link]
                        # stubbify
                        self.agents[name1]['iface'][site1]['bond'] = ''.join([name, self.bond_sep, site])
                        self.agents[name]['iface'][site]['bond'] = ''.join([name1, self.bond_sep, site1])
                    else:
                        pending[link] = [(name, site)]
                        continue
                elif self.bond_sep in link:
                    # this occurs when we created the molecule with an already 'stubbified' agent dictionary
                    degree += 1
                    name1, site1 = link.split(self.bond_sep)
                    complement = self.agents[name1]['iface'][site1]['bond']
                    if complement not in confirmed_stubs:
                        confirmed_stubs += [link]
                        continue
                else:
                    # bond state is a ghost, or '_', or '#'
                    self.is_pattern = True
                    continue
                # We get here only once we've seen both ends of a bond -- the partner of the bond of
                # agent 'name' at site 'site' is 'name1' at site 'site1'. Standardize the bond by
                # sorting, then collect it as a unique key in self.bonds.
                b = tuple(sorted([(name1, site1), (name, site)],
                                 key=lambda x: (alphanum_key(x[0]), alphanum_key(x[1]))))
                self.bonds[b] = 1  # just an indicator; we are collecting unique keys (bonds)
            self.agents[name]['info']['degree'] = degree

    def _stubbify_rebuild(self, id_shift=0, remapping=None):
        """
        Replaces numeric bond labels with unique bond stubs much like _stubbify_inplace(), but also
        relabels agent identifiers (via `remapping`, or by a uniform `id_shift` if no remapping is
        given) and builds a fresh agent dictionary, since the agent names change. Absent a remapping
        and an id_shift, it amounts to _stubbify_inplace() -- prefer that in this case, especially
        for very large complexes, since it avoids rebuilding the whole agent dictionary.

        The label shift is needed when connecting molecules into a larger molecule: the fusion
        requires shifting the agent identifiers of one of the molecules by the number of agents
        contained in the other.
        """
        self.bonds = {}
        pending = {}
        confirmed_stubs = []
        new_agents = {}
        if remapping is None:
            remapping = {entry['info']['id']: str(int(entry['info']['id']) + id_shift)
                         for entry in self.agents.values()}

        for name in self.agents:
            entry = self.agents[name]
            agent_type = entry['info']['type']
            new_id = remapping[entry['info']['id']]
            new_name = add_identifier(agent_type, new_id)
            new_agents[new_name] = {
                'iface': {},
                'info': {'id': new_id, 'type': agent_type, 'sID': entry['info']['sID'], 'degree': 0},
                'local_view': entry['local_view'],  # shift-independent
            }
            degree = 0
            for site, port in entry['iface'].items():
                # bond may be overwritten below, once its partner is found
                new_agents[new_name]['iface'][site] = {'state': port['state'], 'bond': port['bond']}
                link = port['bond']
                if link == '.':
                    self.n_free_sites += 1
                    continue
                if is_number(link):
                    degree += 1
                    if link in pending:
                        name1, site1 = pending[link]
                        new_id1 = remapping[self.agents[name1]['info']['id']]
                        new_name1 = add_identifier(self.agents[name1]['info']['type'], new_id1)
                        new_agents[new_name1]['iface'][site1]['bond'] = ''.join([new_name, self.bond_sep, site])
                        new_agents[new_name]['iface'][site]['bond'] = ''.join([new_name1, self.bond_sep, site1])
                    else:
                        pending[link] = (name, site)
                        continue
                elif self.bond_sep in link:
                    # this occurs when we created the molecule with an already 'stubbified' agent dictionary
                    degree += 1
                    name1, site1 = link.split(self.bond_sep)
                    complement = self.agents[name1]['iface'][site1]['bond']
                    if complement in confirmed_stubs:
                        new_id1 = remapping[self.agents[name1]['info']['id']]
                        new_name1 = add_identifier(self.agents[name1]['info']['type'], new_id1)
                        new_agents[new_name1]['iface'][site1]['bond'] = ''.join([new_name, self.bond_sep, site])
                        new_agents[new_name]['iface'][site]['bond'] = ''.join([new_name1, self.bond_sep, site1])
                    else:
                        confirmed_stubs += [link]
                        continue
                else:
                    # bond state is a ghost, or '_', or '#'
                    self.is_pattern = True
                    continue
                n1 = self.agents[name1]['info']['type']
                (t1, l1, s1), (t2, l2, s2) = sorted([(n1, int(new_id1), site1), (agent_type, int(new_id), site)])
                b = (add_identifier(t1, str(l1)), s1), (add_identifier(t2, str(l2)), s2)
                # collect unique bonds
                self.bonds[b] = 1  # just an indicator; we are collecting unique keys (bonds)
            new_agents[new_name]['info']['degree'] = degree
        self.agents = new_agents

    def get_local_views(self):
        """
        Obtain the lexically ordered local view at each agent.
        """
        # get the last used system-wide index of local views encountered thus far
        if not self.local_view_index:
            running_id = 0
        else:
            running_id = self.local_view_index[next(reversed(self.local_view_index))]

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
            if l_view not in self.local_view_index:
                running_id += 1
                self.local_view_index[l_view] = running_id

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
        if not self.local_view_index:
            return ''
        # get the local view with the smallest index in the _system_ (!)
        # we could also sort lexicographically; but this makes it compatible with
        # subclassing in KappaMolecule.
        _, mlv = min([(self.local_view_index[lv], lv) for lv in self.local_views])
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
            canonic.append(str(self.local_view_index[self.agents[n]['local_view']]))
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
            adjacency = [iface[s1]['bond'].split(self.bond_sep)[0] for s1 in iface if iface[s1]['bond']
                         not in {'.', '#', '_'}]
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
        # wrapper to allow for subclass overriding
        self._remap_ids(remapping)

    def _remap_ids(self, remapping):
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
            self._setup_navigation()

        if self.canon and not self.is_pattern:
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
        Prints a summary of the molecule at various levels of detail.
        """
        dashes = ''.ljust(70, '-')
        info = '\n' + f'{dashes}\n'
        info += self._signature_warning()
        info += f'{self._count_prefix()}{self.size} agents, {len(self.bonds)} bonds, and {self.n_free_sites} free sites\n'
        info += f'composition: {self.sum_formula}\n'
        if self.is_pattern:
            info += 'expression is a pattern !\n'
        info += self.show(internal=internal, wrap=200) + '\n'
        info += f'{dashes}\n'
        if show_bonds:
            info += self.report_bonds(db_level=db_level, pp_width=3)
        info += self._reactivity_report(reactivity, pp_width)
        return info

    def _signature_warning(self):
        """Hook for KappaMolecule: warns when no signature is available for reactivity bookkeeping."""
        return ''

    def _count_prefix(self):
        """Hook for KappaMolecule: prefixes the agent/bond/free-site count line with the instance count."""
        return ''

    def _reactivity_report(self, reactivity, pp_width):
        """Hook for KappaMolecule: appends a reaction-propensity report when requested."""
        return ''

    def report_bonds(self, db_level=0, pp_width=40):
        """
        Lists the bonds of the molecule (verbose, unaware of any signature).
        Overridden by KappaMolecule to additionally report by bond/site type.
        """
        return self._bond_listing(pp_width) if db_level == 2 else ''

    def _bond_listing(self, pp_width):
        info = '\n'
        s = 'list of bonds (random order):'
        info += f'{s:>{pp_width}}\n'
        for (a1, s1), (a2, s2) in self.bonds:
            info += f'{"":>{pp_width}} {a1}{s1}<->{a2}{s2}\n'
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


class KappaMolecule(KappaExpression):
    """
    A KappaExpression plus the bookkeeping needed in the context of a Mixture: instance count,
    a signature-dependent index of free sites and bonds by type (self.free_site_list,
    self.bond_list, and their sibling counters), and the intra-molecular reaction propensities
    derived from them (self.binding, self.unbinding). See KappaExpression for the underlying
    site-graph representation (self.agents, self.adjacency, self.bonds, ...).
    """

    def __init__(self,
                 agents=None,
                 count=0,
                 id_shift=0,
                 signature=None,
                 system=None,
                 local_view_index=None,
                 has_local_views=False,
                 nav=True,
                 canon=True,
                 init=True):

        super().__init__(agents=agents, id_shift=id_shift, nav=nav, canon=canon, init=False)

        self.count = count

        self.system = system
        self.signature = signature

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

        self.has_local_views = has_local_views
        # Local view indices of the mixture in the context of which molecules are canonicalized.
        # Keep in mind that changes in self.local_view_index are implicitly returned to the caller,
        # since a mutable dict is shared, not copied.
        self.local_view_index = {} if local_view_index is None else local_view_index

        # override so we don't have to set signature and local_view_index individually
        if self.system:
            # signature is used for computing internal reaction propensities
            self.signature = self.system.signature
            self.canon = self.system.canonicalize
            self.nav = False
            if self.system.mixture:
                self.local_view_index = self.system.mixture.local_views

        # if data are empty, this generates an empty KappaMolecule
        if not self.agents:
            return

        if init:
            self.initialize_expression(canon=False, nav=self.nav)
            self.initialize_molecule()

    def initialize_molecule(self):
        if self.system:
            if self.size >= self.system.size_threshold:
                # self.canon = False
                self.make_adjacency_lists()
                self.canonical = self.system.served_name
                self.system.served_name += 1

        if self.canon and not self.is_pattern:
            if not self.adjacency:
                self.make_adjacency_lists()
            # requires adjacency list
            if not self.has_local_views:
                self.get_local_views()
            self.make_local_view_lists()
            self.canonical = self.canonicalize()

        # calculate reaction propensities
        if self.signature:
            self._update_reactive_statistics()
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
            self._setup_navigation()

        # calculate reaction propensities
        if self.signature:
            self._update_reactive_statistics()
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

    def _update_reactive_statistics(self):
        """
        (Re)computes signature-dependent bookkeeping used for reaction-propensity calculations:
        per-type free-site and bond lists/counts, and the number of intra-agent binding
        opportunities excluded from the combinatorial bond-formation counts. Assumes bonds are
        already stubbified (self.bonds and self.agents[...]['iface'][...]['bond'] are current).
        """
        self.clear_type_lists()

        for b in self.bonds:
            bt = bond2type(b)
            self.bond_list[bt].append(b)
            self.bond_list_idx[bt][b] = self.bond_type[bt]  # (ab)used as a counter
            self.bond_type[bt] += 1

        for name, entry in self.agents.items():
            agent_type = entry['info']['type']
            # the local free-site types of this agent, used below to detect self-binding opportunities
            agent_free_site_types = set()
            for site, port in entry['iface'].items():
                if port['bond'] == '.':
                    st = f'{agent_type}.{site}'
                    agent_free_site_types.add(st)
                    self.free_site_list[st].append((name, site))
                    self.free_site_list_idx[st][(name, site)] = self.free_site[st]  # (ab)used as a counter
                    self.free_site[st] += 1
            # Accumulate the intra-agent binding opportunities within the whole molecule.
            # This will be used to correct the intra-molecular bond formation propensity
            # (by removing agent self-binding).
            for bt in self.signature.bond_types:
                st1, st2 = bt
                # The case st1 == st2 cannot arise within a single Kappa agent
                if st1 != st2 and st1 in agent_free_site_types and st2 in agent_free_site_types:
                    self.agent_self_binding[bt] += 1

    def remap_ids(self, remapping):
        """
        (Re)assigns agent labels (identifiers) using the map 'remapping'.
        Overrides parent class method
        """
        self._remap_ids(remapping)

        if self.signature:
            self._update_reactive_statistics()

    def _signature_warning(self):
        return '' if self.signature else 'Warning: no signature. Counts presume a signature.\n'

    def _count_prefix(self):
        return f'[count: {self.count}] '

    def _reactivity_report(self, reactivity, pp_width):
        return self.report_internal_reaction_propensities(pp_width=pp_width) if reactivity else ''

    def report_bonds(self, db_level=0, pp_width=40):
        """
        Lists the bond and free-site counts by type (requires a signature); falls back to
        the raw bond listing (see KappaExpression.report_bonds) when db_level == 2.
        """
        if self.signature:
            s = f'{len(self.bond_type)} bond types:'
            info = f'{s:>{pp_width}}\n'
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
            info += self._bond_listing(pp_width)
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


# -------------------------------------------------------------------------------------------


if __name__ == '__main__':
    import kasnap as snap
    import kamatch

    # parser
    kappa = Kappa()
    # a simple Kappa string
    s1 = ' A(o[1], p[2] t{p}[3]), B(x[1] y[2] z[.]), C(w[3], y[z.B])'
    agents = kappa.parse(s1)
    c = KappaExpression(agents)
    out = c.kappa_expression()
    print(f"expression:\n{out}")
    out = c.show(internal=True)
    print(f"internal representation:\n{out}")

    print("--------------")

    s2 = " x222667:P(a1[.] a2[.] a3[.] d[1]), x251065:P(a1[.] a2[.] a3[.] d[1])"
    agents = kappa.parse(s2)
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
    #     molecule_ = kappa_to_molecule(kappa.decode(canonical, snapshot.local_views))
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
    x1 = kappa_to_molecule(data)
    print(x1.kappa_expression())
    print(x1.canonical)
    out = kappa.decode(x1.canonical, x1.local_view_index)
    print(out)
    print(f'decoded is isomorphic to original: {SGM.isomorphic(kappa_to_molecule(out), x1)}')
    expression = canonical_to_expression(x1.canonical, x1.local_view_index)
    print(expression)
    print(f'decoded is isomorphic to original: {SGM.isomorphic(expression, x1)}')
