#!/usr/bin/env python
"""
buildwcs.py
Routines to create YAML WCS files to serve as starting points for WCSFit.
"""

import yaml
import sys
import copy
import numpy as np
from gmbpy.decam import ccdnums

def newRings(device, filename='astrorings2.yaml',coeff=1.):
    """
    Return a dict that will serialize a tree ring template.
    Initialize scaling to the value of coeff
    """
    out = {}
    out['Type'] = 'Template'
    out['HasSplit'] = True
    out['XSplit'] = 1024
    out['Parameter'] = coeff
    out['Filename'] = filename
    out['LowTable'] = device + '_low'
    out['HighTable'] = device + '_high'
    return out

def newXEdge(xStart, xEnd, dx):
    """
    Return a dict that serializes a piecewise correction over the given
    range of x.  Put first nonzero point at xStart, last at xEnd, steps
    of size >=dx.
    Two nodes will be added outside this range which must be zero.
    Initialize all nodes to zero.
    """
    nSteps = int(np.floor((xEnd-xStart) / dx))
    argStep = (xEnd - xStart) / nSteps
    argStart = xStart - argStep
    values = [0.,0.]
    for i in range(nSteps):
        values.append(0.)
    out = {}
    out['Type']='Piecewise'
    out['Axis'] = 'X'
    out['ArgStart'] = argStart
    out['ArgStep'] = argStep
    out['Parameters'] = values
    return out

def newComposite(elements):
    """
    Return a dict that serializes a Composite polymap as the sequence
    of elements in the elements list
    """
    out = {}
    out['Type'] = 'Composite'
    out['Elements'] = elements
    return out

def instruments(names, sharedName, polyDict):
    """
    Make dict of WCS elements specifying a set of instruments
    with given names that all have same tree rings and piecewise edges
    using the sharename.
    Get the starting polynomial values for each device from the
    deciphered YAML stored in the dictionary polyDict.
    """

    leftRange = (30,160)
    rightRange = (1888, 2018)
    xStep = 10.

    skip = ['N30','S30']
    all = {}
    for detpos in ccdnums.keys():
        # First, for each device have a common left/right glowing edge
        if detpos in skip:
            continue
        if detpos=='S31':
            xmin = 50
        elif detpos=='N8':
            xmin = 40
        else:
            xmin = leftRange[0]
        leftName = sharedName + '/' + detpos + '/left'
        all[leftName] = newXEdge(xmin, leftRange[1], xStep)

        if detpos in ['S29','S25']:
            xmax = 2008
        elif detpos in ['N15','N25']:
            xmax = 1988
        else:
            xmax = rightRange[1]
        rightName = sharedName + '/' + detpos + '/right'
        all[rightName] = newXEdge(rightRange[0], xmax, xStep)

        # Make a shared treering template for this device
        treeName = sharedName + '/' + detpos + '/rings'
        all[treeName] = newRings(detpos)

        # Every output instrument gets its own polynomial solution
        for inst in names:
            devName = inst + '/' + detpos
            polyName = devName + '/poly'
            elements = [leftName, rightName, treeName, polyName]
            all[devName] = newComposite(elements)
            all[polyName] = copy.deepcopy(polyDict[sharedName + '/' + detpos])

    return all
    
if __name__=='__main__':
    if len(sys.argv) < 4:
        print "Usage: buildwcs.py <inYAML> <inName> <outname1> [outname2] ..."
        print "stdout is a YAML WCS specification for instruments outname*"
        print "inYAML is file holding the initial polynomial terms"
        print "inName is the instrument name for polynomials, and the"
        print "  shared name for the edges and tree rings"
        sys.exit(1)

    names = sys.argv[3:]
    sharedName = sys.argv[2]
    f = open(sys.argv[1])
    inPoly = yaml.load(f)
    f.close()
    
    all = instruments(names, sharedName, inPoly)
    yaml.dump(all, stream=sys.stdout)

    sys.exit(0)
    
        
