#!/usr/bin/env python
"""
Python script for doing photometric fits given a set of catalogs
from starflat sequences.
This version uses existing color catalogs and does only requested bands
"""

steps = ("begin", "config", "fof", "final", "draw","end")
pieces = False # True to fit piecewise edge functions, False for known templates
fixEdges = True # Set True to hold coefficients of edge templates to 1.
usage = "Script to run photometric fitting process\n"\
        "usage: bandPhoto.py <run name> <bands> [beginStep] [endStep]\n"\
        "   <run name> is base name for files in this run. <run name>.yaml\n"\
        "      file should give configuration.\n"\
        "   <bands> is comma-separated list of bands to process.\n"\
        "   [beginStep]  is name of first step in process to execute.Defaults to 'begin'\n"\
        "   [endStep]    is name of last step in process to execute. Defaults to 'end'\n"\
        "Here are the steps: " + str(steps)
import subprocess
import os
import sys

configDir= "/data3/DECAM/PHOTCAL/STARFLATS"
defaultIn = os.path.join(configDir,'config.yaml')
astref = "/data3/DECAM/PHOTCAL/STARFLATS/astref.yaml"
fields = "/data3/DECAM/PHOTCAL/STARFLATS/fields.yaml"
colors = "/data3/DECAM/PHOTCAL/STARFLATS/colors.yaml"
atts = "/data3/DECAM/PHOTCAL/STARFLATS/attributes.yaml"


if len(sys.argv) <3:
    print usage
    sys.exit(1)

run = sys.argv[1]
configIn = run + '.yaml'
bands = sys.argv[2].split(',')

configFile = run + ".config"
fofFile    = run + ".fof"

# Can optionally control which steps of the process are done
def stepNumber(s):
    for i,v in enumerate(steps):
        if s.lower()==v.lower():
            return i
    print "Unknown configuration step:",s
    sys.exit(1)

    
if len(sys.argv)>3:
    startStep = stepNumber(sys.argv[3])
else:
    startStep = stepNumber("begin")
if len(sys.argv)>4:
    endStep = stepNumber(sys.argv[4])
else:
    endStep = stepNumber("end")
    
    
if startStep <= stepNumber("config") <= endStep:
    # Collect all configuration information
    cmd = ['configure.py']
    cmd.append(fields)
    cmd.append(configIn)
    cmd.append(atts)
    cmd.append(astref)
    cmd.append(colors)
    cmd.append(configFile)
    subprocess.call(cmd)

if startStep <= stepNumber("fof") <= endStep:
    # Match objects
    cmd = ["WCSFoF"]
    cmd.append(configFile)
    cmd.append('-matchradius=2')
    cmd.append('-minmatch=3')
    cmd.append('-outName='+fofFile)
    subprocess.call(cmd)


if startStep <= stepNumber("final") <= endStep:
    # Photo fitting
    for b in bands:
        cmd = ['PhotoFit']
        cmd.append(fofFile)
        cmd.append(os.path.join(configDir,'photofit.params'))
        cmd.append('-useInstruments='+b)
        cmd.append('-outCatalog={:s}.{:s}.photo.fits'.format(run,b))
        cmd.append('-outPhotFile={:s}.{:s}.photo'.format(run,b))
        cmd.append('-outPriorFile={:s}.{:s}.prior.out'.format(run,b))
        cmd.append('-deviceModel')
        if pieces:
            cmd.append('Poly 4 + Color Poly 1 + ' \
                       'Template photorings3.yaml DEVICE + ' \
                       'Piecewise X 8 153 8 + ' \
                       'Piecewise X 1897 2048 8')
        else:
            cmd.append('Poly 4 + Color Poly 1 + ' \
                       'Template photorings3.yaml DEVICE + ' \
                       'Template edges.yaml INST/DEVICE_left + ' \
                       'Template edges.yaml INST/DEVICE_right')
            if fixEdges:
                cmd.append('-fixMaps='+b+'/[NS].*_[34]')
        cmd.append('-exposureModel')
        cmd.append('Constant + Color Constant')
        ###cmd.append('Poly 1 + Color Constant')
        ###cmd.append('Poly 2 + Color Constant')
        cmd.append('-colorExposures=gi0640,gi0730,gi1327,gi1900,gi2040')
               
        outfile = open('{:s}.{:s}.photofit.out'.format(run,b), 'w')

        subprocess.call(cmd,stdout=outfile)
        outfile.close()

if startStep <= stepNumber("draw") <= endStep:
    for b in bands:
        # Draw the flat with rings and edges stripped out
        cmd = ('python ' + os.path.join(configDir,'stripRings.py')).split()
        fin = open('{:s}.{:s}.photo'.format(run,b))
        stripped = 'tmp.photo'
        fout = open(stripped,'w')
        subprocess.call(cmd, stdin=fin, stdout=fout)
        fin.close()
        fout.close()
    
        cmd = 'DrawPhoto 32'.split()
        cmd.append(stripped)
        cmd.append(b)
        cmd.append('{:s}.starflat_{:s}.fits'.format(run,b))
        subprocess.call(cmd)
