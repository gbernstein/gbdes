#!/usr/bin/env python
from __future__ import division,print_function

import sys
from glob import glob
import yaml
from yaml import CLoader as Loader
import re

if len(sys.argv)!=2:
    print("Usage: striprings.py <basename>\n"
              "Will find all files matching <basename>.*.photo and read them as YAML\n"
              "formatted photometric solutions.  It will then find all tree-ring and\n"
              "glowing-edge terms in the solutions and null them out, then *overwrite*\n"
              "the source files with these stripped versions.")
    sys.exit(1)
basename = sys.argv[1]

# Make YAML files stripped of rings and edges
files = glob(basename + ".*.photo")

m1 = re.compile(r'.*rings')
m2 = re.compile(r'.*edge')
for fn in files:
    with open(fn) as f:
        root = yaml.load(f,Loader=Loader)
    # replace rings & edges with Identity
    for k in root:
        if m1.match(k) or m2.match(k):
            root[k] = {'Type':'Identity'}
    root['PhotoMapCollection'] = 'Stripped rings and edges'
    # Write back out to YAML
    print("Writing to",fn)
    with open(fn,'w') as f:
      yaml.dump(root, f)
      

sys.exit(0)
