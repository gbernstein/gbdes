#!/usr/bin/env python
''' Program to extract the exposure-independent components of PixelMaps
from the allsf3 astrometric solutions.

It will also be necessary to map certain filter/epoch combinations into
other epochs in cases where a star flat session was missing filters.
'''
from __future__ import division,print_function
import numpy as np
import yaml
import re
import sys

if len(sys.argv)!=3:
    print("Usage: allsf3.getGuts.py <astro file in> <guts file out>")
    sys.exit(1)

with open(sys.argv[1]) as f:
    root = yaml.load(f, Loader=yaml.CLoader)

ringre = re.compile(r'.*rings')
edgere = re.compile(r'.*edge')
shiftre = re.compile(r'.*ccdshift')
polyre = re.compile(r'[grizY].*poly')
colorre = re.compile(r'[grizY].*color')
instrre = re.compile(r'[grizY].*/[NS]\d+')

for k in list(root.keys()):
    if not (ringre.match(k) or edgere.match(k) or shiftre.match(k) \
      or polyre.match(k) or colorre.match(k) or instrre.match(k)):
      root.pop(k)

root['PixelMapCollection'] = "allsf3 camera solution June 2019"

# Now we need to fix a flaw that there are missing filters in 20180103.
# Force these to use the elements of 20171129 solution (ccdshifts might
# be wrong).
root['z20180103/DEVICE'] = {'Type': 'Composite',
    'Elements': ['z/DEVICE/rings',
                'z/DEVICE/lowedge',
                'z/DEVICE/highedge',
                'z/Y1/DEVICE/poly',
                '20171129/DEVICE/ccdshift',
                'z/color']
                }
root['Y20180103/DEVICE'] = {'Type': 'Composite',
    'Elements': ['Y/DEVICE/rings',
                'Y/DEVICE/lowedge',
                'Y/DEVICE/highedge',
                'Y/Y1/DEVICE/poly',
                '20171129/DEVICE/ccdshift',
                'Y/color']
                }

with open(sys.argv[2],'w') as f:
    yaml.dump(root, f, Dumper=yaml.CDumper)

sys.exit(0)
