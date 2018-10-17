from setuptools import setup, find_packages, Extension
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
                [ os.path.join('AdmixtureBayes','src', '_covariance_matrix.pyx')],
                include_dirs=[np.get_include()],
                extra_compile_args = compile_args)
            ]
            


setup(
        name='AdmixtureBayes',
        version='0.1',
        license='Creative Commons Attribution-Noncommercial-Share Alike license',
        long_description=open('README.md').read(),
        ext_modules = cythonize(extensions), 
        zip_safe=False, 
        
        install_requires=[
            "numpy>=1.11.0",
            "scipy>=0.17.0", 
            "matplotlib", 
            "pandas", 
            "pathos", 
            "pygraphviz"
        ],
        
        package_dir = {'': 'AdmixtureBayes'}, 
        packages=['src'],
        
        entry_points={'console_scripts': [
            'AdmixtureBayes = src.__main__:main'
        ]},
        )
        

