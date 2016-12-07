# setup.py to cythonize metaNW_globalPath.pyx
# executed with python3 /Users/evanbiederstedt/Downloads/setup.py build_ext --inplace:wq

from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'meta NW, global path',
  ext_modules = cythonize("/Users/evanbiederstedt/Downloads/metaNW_globalPath.pyx"),
)
