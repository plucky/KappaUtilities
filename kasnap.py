# Walter Fontana, 2022
"""
This module manages the 'mixture' of KappaMolecules.
"""
# population of Kappa complexes

import re
import os
import json

import kamol


def load_and_unpack(kappa_file, system=None, local_views={}, signature=None, canon=True):
    """
    Load a Kappa snapshot file.
    Format:

    // Snapshot [Event: 662579]
    // "uuid" : "452918066"
    %def: "T0" "1.135"

    %init: 1 /*22 agents*/ A(l[1] r[.] p[2]), A(l[3] r[1] p[4]),
      A(l[.] r[3] p[5]), P(a1[5] a2[.] a3[.] d[6]), P(a1[.] a2[7] a3[.] d[6]),
      A(l[.] r[8] p[7]), A(l[8] r[9] p[10]), A(l[9] r[.] p[11]),
      P(a1[11] a2[12] a3[.] d[.]), A(l[.] r[13] p[12]), A(l[13] r[.] p[14]),
      P(a1[10] a2[14] a3[.] d[15]), P(a1[.] a2[.] a3[4] d[15]),
      P(a1[2] a2[16] a3[.] d[17]), A(l[.] r[.] p[16]),
      P(a1[.] a2[.] a3[18] d[17]), A(l[19] r[.] p[18]), A(l[.] r[19] p[20]),
      P(a1[21] a2[.] a3[20] d[22]), A(l[.] r[.] p[21]),
      P(a1[.] a2[23] a3[.] d[22]), A(l[.] r[.] p[23])
    %init: 1 /*14 agents*/ A(l[.] r[1] p[.]), A(l[1] r[.] p[2]),
      etc.
    """
    complexes = []
    rg_state = None
    lv = local_views  # note: if there is a system.mixture, KappaMolecule will use the mixture's local views
    kappa = kamol.Kappa()  # the parser
    if not os.path.isfile(kappa_file):
        raise Exception("Cannot find snapshot file %s" % kappa_file)
    else:
        with open(kappa_file, "r") as data:
            event = int(data.readline().split('Event:')[1][:-2].strip())
            uuid = data.readline().split('"uuid" : ')[1][1:-2]
            for i in range(0, 4):
                peek = data.readline()
                if "RG" in peek:
                    # RNG state
                    rg_state = json.loads(peek.split('// RG state : ')[1][:-1])
                elif "T0" in peek:
                    t = peek.split('T0')[1][:-2]
                    time = float(re.sub(r'"', ' ', t).strip())
                elif "LV" in peek:
                    # local views included
                    lv = json.loads(peek.split('// LV : ')[1][:-1])
                else:
                    break
            # data.readline()

            while True:
                entry = next_complex_from_file(data)
                if not entry:
                    break
                # parse the entry
                match = re.findall(r'%init: (.*?) \/\*(.*?) agents\*\/ (.*?)$', entry)[0]
                # build the internal representation
                komplex = kamol.KappaMolecule(kappa.parser(match[2].strip()),
                                              count=int(match[0]),
                                              system=system,
                                              sig=signature,
                                              canon=canon,
                                              s_views=lv)  # local_views 'lv' will be updated
                complexes.append(komplex)
    del kappa
    return event, uuid, rg_state, time, complexes, lv


def next_complex_from_file(data):
    """
    Read an %init entry over multiple lines from a snapshot file.
    """
    entry = ''
    first = True
    while True:
        last_pos = data.tell()
        line = data.readline()
        if not line:
            return entry
        if '%init' in line:
            if 'agents*/' in line:
                if not first:
                    data.seek(last_pos)
                    break
                else:
                    first = False
            else:
                return ''
        entry = entry + line[:-1]  # remove \n
    return entry


def gather_local_views_of_snapshot(complexes):
    """
    Collects all local views in a snapshot.
    """
    system_views = {}
    running_id = 0
    for mol in complexes:
        for l_view in mol.local_views:
            if l_view not in system_views:
                running_id += 1
                system_views[l_view] = running_id
    return system_views


def consolidate(complexes):
    """
    Determines equality between complexes and consolidates.
    Requires canonical forms.
    """
    to_delete = []
    n_complexes = len(complexes)
    for i in range(0, n_complexes):
        m1 = complexes[i]
        for j in range(i+1, n_complexes):
            m2 = complexes[j]
            if m2.count == 0:
                continue
            if m1.size == m2.size:
                if m1.sum_formula == m2.sum_formula:
                    if m1.canonical == m2.canonical:
                        # the molecular species exists
                        m1.count += m2.count
                        m2.count = 0
                        to_delete.append(m2)

    for m in to_delete:
        complexes.remove(m)


class SnapShot:
    """
    Reads a snapshot from a .ka file and stores complexes in the list self.complexes of KappaMolecule objects.
    Alternatively import an already-made list of KappaMolecules (complexes).
    """

    def __init__(self, file=None, complexes=[], system=None, signature=None, canon=True):
        self.file = file
        self.origin_uuid = ''
        self.rg_state = None
        self.time = 0
        self.event = 0
        self.number_of_species = 0
        self.complexes = []  # list of 'KappaMolecule's
        # This is a mapping from local views to a running number to determine the 'smallest' local view,
        # according to this mapping, contained in a molecule. This global ordering is needed for the
        # canonicalization of expressions.
        self.local_views = {}

        if file:
            if self.file.endswith(".ka"):
                value = load_and_unpack(file,
                                        system=system,
                                        local_views=self.local_views,
                                        signature=signature,
                                        canon=canon)
                self.event, self.origin_uuid, self.rg_state, self.time, self.complexes, self.local_views = value
            else:
                raise Exception("Unknown file extension %s" % self.file)
        elif complexes:
            self.complexes = complexes
            # get the list of local views in the mixture
            self.local_views = gather_local_views_of_snapshot(self.complexes)

        self.number_of_species = len(self.complexes)

    def get_size_distribution(self, dictionary=False):
        """
        Generates the size distribution present in the mixture.

        Returns a sorted list of pairs (length, count)
                 eg. [(1, 8583), (2, 2642), (7, 836), (4, 253), (103, 82)]
        """
        length = {}
        for c in self.complexes:
            if c.size in length.keys():
                length[c.size] += c.count
            else:
                length[c.size] = c.count

        if dictionary:
            d = {'size': [], 'count': []}
            for l, c in sorted(length.items(), key=lambda i: i[0], reverse=False):
                d['size'].append(l)
                d['count'].append(c)
            return d
        else:
            return [(l, c) for l, c in sorted(length.items(), key=lambda i: i[0], reverse=False)]

    def __str__(self):
        self.snap_report()

    def count_agents_and_molecules(self):
        """
        Counts the number of molecules in the mixture and the total number of agents of each type.
        """
        agents = {}
        molecules = 0
        for m in self.complexes:
            molecules += m.count
            for a in m.composition:
                if a in agents:
                    agents[a] += m.composition[a] * m.count
                else:
                    agents[a] = m.composition[a] * m.count
        return agents, molecules

    def snap_report(self, pp_width=40):
        """
        Summarizes the mixture.
        """
        info = f"\n{'MIXTURE '.ljust(70, '-')}\n\n"
        info += f'{"initial mixture file":>20}: {self.file}\n'
        info += f'{"molecular species":>20}: {self.number_of_species}\n'
        info += f'{"agents, molecules":>20}: {self.count_agents_and_molecules()}\n'

        info += f'{"size distribution":>20}:\n'
        size_dist = self.get_size_distribution()
        d1 = f'{size_dist[:4]}'
        d2 = f'{size_dist[-4:]}'
        info += f'{d1[1:-1]} ... {d2[1:-1]}'
        info += '\n'
        return info
        
# -------------------------------------------------------------------------------------------


if __name__ == '__main__':

    snap = SnapShot('TestData/snap__1773.ka')
    print(f'number of species: {snap.number_of_species}')
    print(f'number of local views at system level: {len(snap.local_views)}')
    print("size distribution")
    print(snap.get_size_distribution())
    print(f'agents and molecules: {snap.count_agents_and_molecules()}')
    print("All complexes share the system-wide local views:")
    claim = True
    first = True
    for m in snap.complexes:
        if first:
            first = False
        else:
            if m.system_views != m_previous.system_views:
                print("Claim is False.")
                break
        m_previous = m
    print("Claim is True.")

