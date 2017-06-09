from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os

# extra compiler flags to make it highly optimized and fast. `-march=native` is
# a gcc-specific flag, I think.
compile_args = ['-march=native', '-O3']

extensions = [
            Extension(
                # oddly, the module name (here, _interp) has to be the same as
                # the basename (ie prefix) of the full filename.
                '_covariance_matrix',
                ['_covariance_matrix.pyx'],
                include_dirs=[np.get_include()],
                extra_compile_args = compile_args)
            ]

setup(
        ext_modules = cythonize(extensions)
        )
