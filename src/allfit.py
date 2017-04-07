#!/usr/bin/env python
"""
Python script for doing photometric fits given a set of catalogs
from starflat sequences.
"""
import subprocess
import os
import os.path
import sys

steps = ("configure", "dcr", "recenter", "match", "photo", "color","astro")

parser = argparse.ArgumentParser(description='Derive photometric and astrometric solutions from a set of detection catalogs')
parser.add_argument("basename", help='Root name for solution-specific files', type=str)
parser.add_argument('--setup_dir', help='Directory holding generic setup files', type=str)
parser.add_argument('--first',help='First step of process to execute',type=str,choices=steps,
                        default=steps[0])
parser.add_argument('--last',help='Last step of process to execute',type=str, choices=steps,
                        default=step[-1])
parser.add_argument('--bands',help='Comma-separated list of bands for photometry', type=str,
                        default='g,r,i,z,Y')

args = parser.parse_args()

# First a little routine that keeps track of whether we are in the range of first--last steps
stepsActive = False   # True if we are to currently executing steps
stepsDone = False   # True if we are done with all requested steps
def doStep(stepname):
    if stepsDone:
        sys.exit(0)  # past last step
    if stepname==args.first:
        stepsActive = True
    if stepname==args.last:
        stepsDone = False
    return stepsActive

# Types of files that will be expected to exist locally:
specFile = "{:s}.setup".format(parser.basename)
priorFile = "{:s}.photo.prior".format(parser.basename)

# Parameter and configuration files
stdSpecs = os.path.join(args.setup_dir,'std.specs')
fofParams = os.path.join(args.setup_dir,'wcsfof.params')
photoParams = os.path.join(args.setup_dir,'photofit.params')
magcolorParams = os.path.join(args.setup_dir,'magcolor.params')
astroParams = os.path.join(args.setup_dir,'wcsfit.params')
# Model specifications
photoInput = os.path.join(args.setup_dir,'nocolor.input.photo')
colorInput = os.path.join(args.setup_dir,'color.input.photo')
astroInput = os.path.join(args.setup_dir,'input.astro')
# Pre-solved model components
oldPhoto = os.path.join(args.setup_dir,'allsf.guts.photo')
oldAstro = os.path.join(args.setup_dir,'allsf.guts.astro')

# Files we will create
configFile = args.basename + '.config'
configLog = args.basename + '.config.log'
dcrFile = args.basename + '.dcr.astro'
fofFile = args.basename + '.fof'
fofLog = args.basename + '.fof.log'
photoFile = args.basename + '.{:s}.photo'
photoCat = args.basename + '.{:s}.photo.fits'
photoLog = args.basename + '.{:s}.photofit.log'
priorLog = args.basename + '.{:s}.prior.out'
colorFile = args.basename + '.gi.cat'
magFile = args.basename + '.mags'
magLog = args.basename + '.magcolor.log'
colorFofFile = args.basename + '.color.fof'
astroFile = args.basename + '.astro'
astroLog = args.basename + '.wcsfit.log'
astroCat = args.basename + '.astro.fits'


bands = set(args.bands.split(','))

if doStep('config'):    
    # Collect all configuration information
    cmd = ['configure.py']
    cmd.append(specfile)
    if os.path.isfile(stdSpecs):
        # More attributes, etc., if given:
        cmd.append(stdSpecs)
    cmd.append(configFile)
    with open(configLog) as log:
        subprocess.check_call(cmd, stdout=log, stderr=subprocess.STDOUT)

if doStep('dcr'):    
    # Make DCR file
    cmd = ['decamDCR.py']
    cmd.append(configFile)
    cmd.append(dcrFile)
    subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        
if doStep('recenter'):    
    # Recenter exposures' RA/Dec
    cmd = ['recenter.py']
    cmd.append(configFile)
    subprocess.check_call(cmd, stdout=sys.stdout, stderr=sys.stderr)
        
if doStep("match"):
    # Match objects
    cmd = ["WCSFoF"]
    cmd.append(configFile)
    cmd.append(fofParams)
    cmd.append('-outName='+fofFile)
    with open(fofLog) as log:
        subprocess.check_call(cmd, stdout=log, stderr=subprocess.STDOUT)

if doStep("photo"):
    # Photo fitting, without color - do g and i
    # ???
    for b in bands & set(['g','i']):
        cmd = ['PhotoFit']
        cmd.append(photoParams)
        cmd.append('-useInstruments='+b+'.*')  # ?? more complex instruments ??
        cmd.append('-outCatalog='+photoCat.format(b))
        cmd.append('-outPriorFile='+priorLog.format(b))
        cmd.append('-outPhotFile='+photoCat.format(b))
        cmd.append('-colorExposures=""')
        if os.path.isfile(oldPhoto):
            cmd.append('-inputMaps='+photoInput+','+oldPhoto)
        else:
            cmd.append('-inputMaps='+photoInput)
        with open(photoLog.format(b)) as log:
            subprocess.check_call(cmd, stdout=log, stderr=subprocess.STDOUT)

if doStep("color"):
    # MagColor, refit, MagColor again...
    pass

# ?? Phot2DESDM

if doStep('astro')
    pass

sys.exit(0)
