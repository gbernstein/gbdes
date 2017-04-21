#!/usr/bin/env python
"""
Python script for doing photometric fits given a set of catalogs
from starflat sequences.
"""
import subprocess
import os
import os.path
import sys
import argparse

steps = ("config", "dcr", "recenter", "match", "prior", "photo", "astro", "draw")

parser = argparse.ArgumentParser(description='Derive photometric and astrometric solutions from a set of detection catalogs')
parser.add_argument("basename", help='Root name for solution-specific files', type=str)
parser.add_argument('--setup_dir', help='Directory holding generic setup files', type=str)
parser.add_argument('--color', help='Produce colors and color terms',
                    type=bool, default=True)
parser.add_argument('--first',help='First step of process to execute',type=str,choices=steps,
                        default=steps[0])
parser.add_argument('--last',help='Last step of process to execute',type=str, choices=steps,
                        default=steps[-1])
parser.add_argument('--bands',help='Comma-separated list of bands for photometry', type=str,
                        default='g,r,i,z,Y')
parser.add_argument('--prior',help='Names of photometric prior files, if pre-constructed', type=str)

args = parser.parse_args()

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
astroParams = os.path.join(args.setup_dir,'wcsfit.params')
# Model specifications
photoInput = os.path.join(args.setup_dir,'nocolor.input.photo')
colorInput = os.path.join(args.setup_dir,'color.input.photo')
astroInput = os.path.join(args.setup_dir,'color.input.astro')
astroNocolorInput = os.path.join(args.setup_dir,'nocolor.input.astro')
# Pre-solved model components
oldPhoto = os.path.join(args.setup_dir,'allsf.guts.photo')
oldAstro = os.path.join(args.setup_dir,'allsf.guts.astro')

# Files we will create
dbFile = args.basename + '.db'
configLog = args.basename + '.config.log'
dcrFile = args.basename + '.dcr.astro'
gradientFile = args.basename + '.gradient.photo'
fofFile = args.basename + '.fof'
fofLog = args.basename + '.fof.log'
if args.prior:
    priorFile = args.prior   # Use pre-existing of this name
else:
    priorFile = args.basename + ".{:s}.prior"

photoFile = args.basename + '.{:s}.photo'
photoCat = args.basename + '.{:s}.photo.fits'
photoLog = args.basename + '.{:s}.photofit.log'
priorLog = args.basename + '.{:s}.prior.out'
colorFile = args.basename + '.color.cat'
magFile = args.basename + '.mags.cat'
magLog = args.basename + '.magcolor.log'
colorFofFile = args.basename + '.color.fof'
astroFile = args.basename + '.astro'
astroLog = args.basename + '.wcsfit.log'
astroCat = args.basename + '.astro.fits'
starflatName = args.basename + '.{:s}.sf'

bands = set(args.bands.split(','))

if doStep('config'):
    print 'Doing config'
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
    # Make DCR file
    print 'Doing dcr'
    cmd = ['decamDCR.py']
    cmd.append(dbFile)
    cmd.append(dcrFile)
    cmd.append(gradientFile)
    subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        
if doStep('recenter'):    
    # Recenter exposures' RA/Dec
    print 'Doing recenter'
    cmd = ['recenter.py']
    cmd.append(dbFile)
    subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        
if doStep("match"):
    # Match objects
    print 'Doing WCSFoF'
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
    if args.color:
        colorBands = ['g','i']
        for b in colorBands:
            if b not in bands:
                print 'Missing band',b,'required for color terms'
                sys.exit(1)

        print 'Doing PhotoFit w/o color'
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
            if os.path.isfile(oldPhoto):
                cmd.append('-inputMaps='+','.join([photoInput,gradientFile,oldPhoto]))
            else:
                cmd.append('-inputMaps='+','.join([photoInput,gradientFile]))
            with open(photoLog.format(b),'w') as log:
                subprocess.check_call(cmd, stdout=log, 
                                      stderr=subprocess.STDOUT)

        print 'Doing first MagColor'
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

        print 'Doing PhotoFit w/ color'
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
            if os.path.isfile(oldPhoto):
                cmd.append('-inputMaps='+','.join([colorInput,gradientFile,oldPhoto]))
            else:
                cmd.append('-inputMaps='+','.join([colorInput,gradientFile]))
            with open(photoLog.format(b),'w') as log:
                subprocess.check_call(cmd, stdout=log, 
                                      stderr=subprocess.STDOUT)

        print 'Doing second MagColor'
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

    # Now make final photofits, with or without color terms
    print 'Doing final PhotoFit'
    for b in bands:
        cmd = ['PhotoFit']
        if args.color:
            cmd.append(colorFofFile)
        else:
            cmd.append(fofFile)
        cmd.append(photoParams)
        cmd.append('-useInstruments='+b+'.*') 
        cmd.append('-priorFiles='+priorFile)
        cmd.append('-outCatalog='+photoCat.format(b))
        cmd.append('-outPriorFile='+priorLog.format(b))
        cmd.append('-outPhotFile='+photoFile.format(b))
        if args.color:
            cmd.append('-colorExposures='+colorFile)
            if os.path.isfile(oldPhoto):
                cmd.append('-inputMaps='+','.join([colorInput,oldPhoto,gradientFile]))
            else:
                cmd.append('-inputMaps='+','.join([colorInput,gradientFile]))
        else:
            cmd.append('-colorExposures=""')
            if os.path.isfile(oldPhoto):
                cmd.append('-inputMaps='+','.join([photoInput,oldPhoto,gradientFile]))
            else:
                cmd.append('-inputMaps='+','.join([photoInput,gradientFile]))

        with open(photoLog.format(b),'w') as log:
            subprocess.check_call(cmd, stdout=log, 
                                  stderr=subprocess.STDOUT)

if doStep('astro'):
    print 'Doing WCSFit'
    cmd = ['WCSFit']
    if args.color:
        cmd.append(colorFofFile)
        cmd.append(astroParams)
        cmd.append('-inputMaps=' + ','.join([oldAstro,astroInput,dcrFile]))
        cmd.append('-colorExposures='+colorFile)
    else:
        cmd.append(fofFile)
        cmd.append('-inputMaps=' + ','.join([oldAstro,astroNocolorInput]))
        cmd.append('-colorExposures=""')
    cmd.append('-outcatalog='+astroCat)
    cmd.append('-outwcs='+astroFile)
    with open(astroLog,'w') as log:
        subprocess.check_call(cmd, stdout=log, 
                              stderr=subprocess.STDOUT)

if doStep('draw'):
    for b in bands:
        print 'Doing Phot2DESDM band',b
        cmd = ['Photo2DESDM']
        cmd.append(photoFile.format(b))
        cmd.append(b) # ??? Here assuming that the instrument name = band
        cmd.append(starflatName.format(b))
        subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)

    print 'Doing DrawAstro'

sys.exit(0)
