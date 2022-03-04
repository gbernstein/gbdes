""" 
DECam-specific information and tasks
"""

# Put utilities into current namespace
from . decamUtils import plotColors,amps,dps,overscanColumns,centralAmpPix,allAmpPix
from . decamUtils import ccdnums, ccdCorners, ccdBounds, ccdCenters, minimalHeader
from . decamUtils import XtalkRemover
from . decamTrig import parallactic, dcr, recenter
from . epochFinder import mjdOfEpoch, EpochFinder

