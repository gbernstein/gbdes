from __future__ import print_function
import sys,os,glob,re
import select

from distutils.core import setup
import distutils

print('Python version = ',sys.version)
py_version = "%d.%d"%sys.version_info[0:2]  # we check things based on the major.minor version.

dependencies = ['numpy', 'astropy']

# Read in the version from gbutil/_version.py
# cf. http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
version_file=os.path.join('gbutil','_version.py')
verstrline = open(version_file, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    gbutil_version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (version_file,))
print('gbutil version is %s'%(gbutil_version))

dist = setup(
        name="gbutil",
        version=gbutil_version,
        author="Gary Bernstein",
        author_email="garyb@PHYSICS.UPENN.EDU",
        description="Python utilities and DECam configuration information",
        license = "BSD License",
        url="https://github.com/gbernstein/gbutil",
        download_url="https://github.com/gbernstein/gbutil/releases/tag/v%s.zip"%gbutil_version,
        packages=['gbutil','gbutil.decam'],
        install_requires=dependencies,
    )

