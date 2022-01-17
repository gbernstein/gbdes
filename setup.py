from pybind11.setup_helpers import Pybind11Extension, build_ext
import os
from glob import glob
from setuptools import setup

yaml_path = os.path.join(os.environ['YAML_DIR'], 'include')
eigen_path = os.environ['EIGEN_DIR']
fftw_path = os.path.join(os.environ['FFTW_DIR'], 'include/')
cfitsio_path = os.path.join(os.environ['CFITSIO_DIR'], 'include/')
include_dirs = [yaml_path, fftw_path, eigen_path, cfitsio_path,
                                    'include',
                                    '../gb_packages/gbfits'] + glob('../gb_packages/*/include')

ext_modules = [
    Pybind11Extension("wcsfit",
                      sources=["pydir/wcsfit.cpp"],
                      define_macros=[('USE_EIGEN', True),
                                     ('USE_YAML', True)],
                      include_dirs=include_dirs,
                      libraries=['gbdes', 'yaml-cpp', 'cfitsio', 'fftw3'],
                      library_dirs=['obj', '../yaml-cpp/lib',
                                    os.path.join(os.environ['CONDA_PREFIX'], 'lib')],
                      extra_compile_args=['-w'],
    )
]

setup(ext_modules=ext_modules)
