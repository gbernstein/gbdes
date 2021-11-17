from pybind11.setup_helpers import Pybind11Extension, build_ext
import os
from glob import glob
from setuptools import setup

yaml_path = os.path.join(os.environ['YAML_DIR'], 'include')
eigen_path = '/software/lsstsw/stack_20210813/stack/miniconda3-py38_4.9.2-0.7.0/Linux64/eigen/3.3.7.lsst2-1-g398bedf/include/'
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
                      libraries=['gbdes', 'yaml-cpp', 'cfitsio', 'fftw3'],#, 'eigen'],
                      library_dirs=['obj', '/home/csaunder/software/yaml-cpp/lib',
                                    os.path.join(os.environ['CONDA_PREFIX'], 'lib')],
                      #library_dirs=['build/temp.linux-x86_64-3.8a/']
                      extra_compile_args=['-w'],#['-Wno-sign-compare', '-Wno-reorder']
    )
]

setup(ext_modules=ext_modules)
