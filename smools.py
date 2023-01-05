# Walter Fontana 2022
"""
This module contains small tools...
"""


def complex2set(cmplx):
    """
    Returns a KappaComplex just as a pair (agent_set, bond_set). Requires sID for each agent.
    """
    agent_set = set()
    bond_set = set()

    if len(cmplx.bonds) == 0:
        for a in cmplx.agents:
            agent_set.add(a['info']['sID'] + a['info']['type'])
    else:
        for (agent1, site1), (agent2, site2) in cmplx.bonds:
            a1 = cmplx.agents[agent1]['info']['sID'] + cmplx.agents[agent1]['info']['type']
            a2 = cmplx.agents[agent2]['info']['sID'] + cmplx.agents[agent2]['info']['type']
            agent_set.add(a1)
            agent_set.add(a2)
            b = sorted([(a1, site1), (a2, site2)], key=lambda i: i[0])
            bond_set.add(tuple(b))

    return agent_set, bond_set


def cut_and_cup(c1, c2):
    """
    returns intersection and union of agent and bond set
    """
    agent_set1, bond_set1 = complex2set(c1)
    agent_set2, bond_set2 = complex2set(c2)

    return agent_set1 & agent_set2, agent_set1 | agent_set2, bond_set1 & bond_set2, bond_set1 | bond_set2


def kappa_expression(agents, bonds, sep='@'):
    """
    Convert KappaMolecule-agents and -bonds into a kappa string.
    Returns a kappa expression.

    agents: dictionary of agents KappaMolecule style
    bonds: set of bonds
    sep: bond separation symbol used by KappaMolecule
    """
    i = 1
    num = {}
    for (p, q) in bonds:
        num[p] = i
        num[q] = i
        i += 1

    s = ''
    for name in agents:
        s += f'{name.split(".")[0]}('
        for site in agents[name]:
            s += f'{site}'
            state = agents[name][site]['state']
            if state != '#':
                s += '{' + f'{state}' + '}'
            link = agents[name][site]['bond']
            if link == '.':
                s += '[.] '
            elif link == '#':
                # s += '[#] '
                s += ''
            elif sep in link:
                ag, ste = link.split(sep)
                s += f'[{num[(ag, ste)]}] '
            else:
                s += f'[{link}] '
        s = s[:-1] + '), '
    return s[:-2]
