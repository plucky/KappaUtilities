# Walter Fontana at 8/5/26

import os
import runpy
import re
import zlib
import base64
from pathlib import Path

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

def pykappa_state_to_snapshot(system_state_file=None, output_directory=None, delta_t=0, with_id=False, sort=True):
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

    def agent_count(line):
        match = re.search(r"/\*(\d+)\s*agents?\*/", line)
        return int(match.group(1))

    system = pykappa.System.load(system_state_file)
    events = sum([r.applied for r in system.tallies.values()])
    uuid = f'"0000"'

    p = Path(system_state_file)

    if not delta_t:
        time = system.time
    else:
        number_str = p.stem.split("_")[-1]
        number = int(number_str)
        # adjust snap.time
        time = number * delta_t

    if with_id:
        lines = system.mixture.kappa_str_with_agent_ids.split('\n')
    else:
        lines = system.mixture.kappa_str.split('\n')
    lines = list(map(add_agent_count, lines))  # add size
    if sort:  # sort on size
        lines = sorted(lines, key=agent_count, reverse=True)

    ka_file = f'{output_directory}/{p.stem}.ka'
    with open(ka_file, "w") as fp:
        prefix = f'// Snapshot [Event: {events}]\n// "uuid" : {uuid}\n%def: "T0" "{time}"\n\n'
        fp.write(prefix)
        for l in lines:
            fp.write(l)
    return ka_file


def binary_decompress(compressed_text):
    # inverse of binary_compress()
    #
    compressed_bytes = base64.b64decode(compressed_text.encode('utf-8'))
    decompressed_bytes = zlib.decompress(compressed_bytes)
    return decompressed_bytes.decode('utf-8')


def binary_compress(text):
    # Convert to bytes
    text_bytes = text.encode('utf-8')
    compressed_bytes = zlib.compress(text_bytes, level=6)
    #  safe ASCII string using Base64
    return base64.b64encode(compressed_bytes).decode('utf-8')


def get_maximer(system):
    return max((component for component in system.mixture), key=lambda c: len(c))
