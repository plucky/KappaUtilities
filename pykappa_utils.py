# Walter Fontana at 8/5/26

import os
import runpy
import re
import pykappa

def start_pykappa_system(system_definition_file=None):
    """
    Start a pykappa system.
    """
    if os.path.isfile(system_definition_file):
        # dynamic execution
        print(f'Using system definition on file ({system_definition_file})')
        func = runpy.run_path(system_definition_file)
        system_ka, p = func["define_system"]()
        system = pykappa.System.from_ka(system_ka, seed=p['seed'])
        return system
    else:
        print(f'system definition file {system_definition_file} not found.')
        return None

def pykappa_state_to_snapshot(system_state=None, temp_directory=None):
    """
    Convert a pykappa system.save() to a KaSim-style snapshot readable by kasnap.py
    """
    # define and start a pykappa system to access the saved system state.
    def add_agent_count(line):
        def replacer(m):
            prefix, num, rest = m.group(1), m.group(2), m.group(3)
            count = rest.count('(')
            return f"{prefix}{num} /*{count} agents*/ {rest}\n"
        return re.sub(r'(%init:\s*)(\d+)\s*(.*)', replacer, line)

    system = pykappa.System.load(system_state)
    events = sum([r['applied'] for r in system.tallies.values()])
    uuid = f'"0000"'
    time = system.time
    temp_file = f'{temp_directory}/temp.ka'

    with open(temp_file, "w") as fp:
        prefix = f'// Snapshot [Event: {events}]\n// "uuid" : {uuid}\n%def: "T0" "{time}"\n\n'
        fp.write(prefix)
        for l in system.mixture.kappa_str.split('\n'):
            fp.write(add_agent_count(l))
    return temp_file
