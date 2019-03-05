#!/usr/bin/env python
"""
Python script for doing photometric fits given a set of catalogs
from starflat sequences.
"""
from __future__ import division,print_function

import subprocess
import os
import os.path
import sys
import argparse
import yaml
from glob import glob

steps = ("catprep","config", "dcr", "recenter", "match", "prior", "photo", "astro", "draw")

parser = argparse.ArgumentParser(description='Derive photometric and astrometric solutions from a set of detection catalogs')
parser.add_argument("basename", help='Root name for solution-specific files', type=str)
parser.add_argument('--setup_dir', help='Directory holding generic setup files', type=str)
parser.add_argument('--no-color', help='Skip colors and color terms',
                    action='store_true')
parser.add_argument('--first',help='First step of process to execute',type=str,choices=steps,
                        default=steps[0])
parser.add_argument('--last',help='Last step of process to execute',type=str, choices=steps,
                        default=steps[-1])
parser.add_argument('--bands',help='Comma-separated list of bands for photometry', type=str,
                        default='g,r,i,z,Y')
parser.add_argument('--prior',help='Names of photometric prior files, if pre-constructed', type=str)
parser.add_argument('--dcr_file',help='Name of pre-existing DCR file; "omit" to omit DCR and gradient corrections', type=str)
parser.add_argument('--skip_recenter',help='Skip exposure recentering',action='store_true')

args = parser.parse_args()

do_color = not args.no_color

# First a little routine that keeps track of whether we are in the range of first--last steps
def doStep(stepname):
    if doStep.stepsDone:
        sys.exit(0)  # past last step
    if stepname==args.first:
        doStep.stepsActive = True
    if stepname==args.last:
        doStep.stepsDone = True
    return doStep.stepsActive
doStep.stepsActive = False   # True if we are to currently executing steps
doStep.stepsDone = False   # True if we are done with all requested steps

# Types of files that will be expected to exist locally:
configFile = "{:s}.config".format(args.basename)

# Parameter and configuration files
stdConfig = os.path.join(args.setup_dir,'std.config')
fofParams = os.path.join(args.setup_dir,'wcsfof.params')
photoParams = os.path.join(args.setup_dir,'photofit.params')
magcolorParams = os.path.join(args.setup_dir,'magcolor.params')
# Model specifications
photoInput = os.path.join(args.setup_dir,'nocolor.input.photo')
colorInput = os.path.join(args.setup_dir,'color.input.photo')
# Pre-solved model components
oldPhoto = os.path.join(args.setup_dir,'allsf.guts.photo')

# Files we will create
dbFile = args.basename + '.db'
configLog = args.basename + '.config.log'

gradientFile = args.basename + '.gradient.photo'
if not args.dcr_file:
    makeDCR = True
    dcrFile = args.basename + '.dcr.astro'
elif args.dcr_file=='omit':
    makeDCR = False
    dcrFile = None
else:
    makeDCR = False
    dcrFile = args.dcr_file
fofFile = args.basename + '.fof'
fofLog = args.basename + '.fof.log'
if args.prior:
    priorFile = args.prior   # Use pre-existing of this name
else:
    priorFile = args.basename + ".{:s}.prior"

photoFile = args.basename + '.{:s}.photo'
photoCat = args.basename + '.{:s}.photo.fits'
photoLog = args.basename + '.{:s}.photofit.log'
photoLog_2 = args.basename + '.{:s}.photofit_2.log'
photoLog_3 = args.basename + '.{:s}.photofit_3.log'
priorLog = args.basename + '.{:s}.prior.out'
colorFile = args.basename + '.color.cat'
magFile = args.basename + '.mags.cat'
magLog = args.basename + '.magcolor.log'
magLog_2 = args.basename + '.magcolor_2.log'
colorFofFile = args.basename + '.color.fof'
starflatName = args.basename + '.{:s}.sf'

bands = set(args.bands.split(','))

if doStep('catprep'):
    print('Doing catalog preparation')
    print('#EXPNUM  MJD_MID    EXPTIME   DEC            HA      BAND APCORR (RMS) FLUXRAD (RMS)')

    # Run catalog preparation on every file
    # matching the globs in the input config
    with open(configFile,'r') as f:
        root = yaml.load(f)
        if 'Files' not in root:
            print('Missing "Files" key in the config file')
            exit(1)
        for d in root['Files']:
            if 'glob' not in d:
                print('Missing "glob" key in FILES spec')
                exit(1)
            for cat in glob(d['glob']):
                subprocess.check_call('frankcat.py '.split() + [cat])

if doStep('config'):
    print('Doing config')
    # Collect all configuration information
    cmd = ['configure.py']
    if os.path.isfile(stdConfig):
        # More attributes, etc., if given:
        cmd.append(stdConfig)
    # Run-specific configs get higher priority by coming later:
    cmd.append(configFile)
    cmd.append(dbFile)
    with open(configLog,'w') as log:
        subprocess.check_call(cmd, stdout=log, stderr=subprocess.STDOUT)

if doStep('dcr'):    
    if makeDCR:
        # Make DCR file
        print('Doing dcr')
        cmd = ['decamDCR.py']
        cmd.append(dbFile)
        cmd.append(dcrFile)
        cmd.append(gradientFile)
        subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        
if doStep('recenter') and not args.skip_recenter:    
    # Recenter exposures' RA/Dec
    print('Doing recenter')
    cmd = ['recenter.py']
    cmd.append(dbFile)
    subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        
if doStep("match"):
    # Match objects
    print('Doing WCSFoF')
    cmd = ["WCSFoF"]
    cmd.append(dbFile)
    cmd.append(fofParams)
    cmd.append('-outName='+fofFile)
    with open(fofLog,'w') as log:
        subprocess.check_call(cmd, stdout=log, stderr=subprocess.STDOUT)
    

if doStep("prior"):
    if not args.prior:
        for b in bands:
            cmd = ['nightlyPrior.py']
            cmd.append(fofFile)
            cmd.append(b)
            cmd.append(priorFile.format(b))
            subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)

if doStep("photo"):
    if do_color:
        colorBands = ['g','i']
        for b in colorBands:
            if b not in bands:
                print('Missing band',b,'required for color terms')
                sys.exit(1)

        print('Doing PhotoFit w/o color')
        for b in colorBands:
            cmd = ['PhotoFit']
            cmd.append(fofFile)
            cmd.append(photoParams)
            cmd.append('-useInstruments='+b+'.*') 
            cmd.append('-priorFiles='+priorFile.format(b))
            cmd.append('-outCatalog='+photoCat.format(b))
            cmd.append('-outPriorFile='+priorLog.format(b))
            cmd.append('-outPhotFile='+photoFile.format(b))
            cmd.append('-colorExposures=""')
            inmaps = [photoInput]
            if os.path.isfile(gradientFile):
                inmaps.append(gradientFile)
            if os.path.isfile(oldPhoto):
                inmaps.append(oldPhoto)
            cmd.append('-inputMaps='+','.join(inmaps))
            with open(photoLog.format(b),'w') as log:
                subprocess.check_call(cmd, stdout=log, 
                                      stderr=subprocess.STDOUT)

        print('Doing first MagColor')
        colorCmd = '-colors='+colorBands[0]+'-'+colorBands[1]+\
                   '@'+colorFile
        cmd = ['MagColor']
        cmd.append(fofFile)
        cmd.append(colorFofFile)
        cmd.append(magFile)
        cmd.append(magcolorParams)
        cmd.append('-bands='+','.join(colorBands))
        cmd.append('-photofiles='+','.join([photoFile.format(b) for b in colorBands]))
        cmd.append(colorCmd)
        with open(magLog.format(b),'w') as log:
            subprocess.check_call(cmd, stdout=log, 
                                  stderr=subprocess.STDOUT)

        print('Doing PhotoFit w/ color')
        for b in colorBands:
            cmd = ['PhotoFit']
            cmd.append(colorFofFile)
            cmd.append(photoParams)
            cmd.append('-useInstruments='+b+'.*') 
            cmd.append('-priorFiles='+priorFile.format(b))
            cmd.append('-outCatalog='+photoCat.format(b))
            cmd.append('-outPriorFile='+priorLog.format(b))
            cmd.append('-outPhotFile='+photoFile.format(b))
            cmd.append('-colorExposures='+colorFile)
            inmaps = [photoInput]
            if os.path.isfile(gradientFile):
                inmaps.append(gradientFile)
            if os.path.isfile(oldPhoto):
                inmaps.append(oldPhoto)
            cmd.append('-inputMaps='+','.join(inmaps))
            with open(photoLog_2.format(b),'w') as log:
                subprocess.check_call(cmd, stdout=log, 
                                      stderr=subprocess.STDOUT)

        print('Doing second MagColor')
        cmd = ['MagColor']
        cmd.append(fofFile)
        cmd.append(colorFofFile)
        cmd.append(magFile)
        cmd.append(magcolorParams)
        cmd.append('-bands='+','.join(colorBands))
        cmd.append('-photofiles='+','.join([photoFile.format(b) for b in colorBands]))
        cmd.append(colorCmd)
        with open(magLog_2.format(b),'w') as log:
            subprocess.check_call(cmd, stdout=log, 
                                  stderr=subprocess.STDOUT)

    # Now make final photofits, with or without color terms
    print('Doing final PhotoFit')
    for b in bands:
        cmd = ['PhotoFit']
        if do_color:
            cmd.append(colorFofFile)
        else:
            cmd.append(fofFile)
        cmd.append(photoParams)
        cmd.append('-useInstruments='+b+'.*') 
        cmd.append('-priorFiles='+priorFile.format(b))
        cmd.append('-outCatalog='+photoCat.format(b))
        cmd.append('-outPriorFile='+priorLog.format(b))
        cmd.append('-outPhotFile='+photoFile.format(b))
        if do_color:
            cmd.append('-colorExposures='+colorFile)
            inmaps = [colorInput]
        else:
            cmd.append('-colorExposures=""')
            inmaps = [photoInput]
        if os.path.isfile(oldPhoto):
            inmaps.append(oldPhoto)
        if os.path.isfile(gradientFile):
            inmaps.append(gradientFile)
        cmd.append('-inputMaps='+','.join(inmaps))
        with open(photoLog_3.format(b),'w') as log:
            subprocess.check_call(cmd, stdout=log, 
                                  stderr=subprocess.STDOUT)


if doStep('draw'):
    for b in bands:
        print('Doing Phot2DESDM for band',b)
        cmd = ['Photo2DESDM']
        cmd.append(photoFile.format(b))
        cmd.append(b) # ??? Here assuming that the instrument name = band
        cmd.append(starflatName.format(b))
        subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)

print("Done")

sys.exit(0)
