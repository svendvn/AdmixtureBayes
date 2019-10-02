from setuptools import setup, find_packages, Extension
import os

try:
    from Cython.Build import cythonize
    import numpy as np
except ImportError:
    use_cython = False
else:
    use_cython = True

if use_cython:
    # extra compiler flags to make it highly optimized and fast. `-march=native` is
    # a gcc-specific flag, I think.
    compile_args = ['-march=native', '-O3']

    extensions = cythonize([
                Extension(
                    # oddly, the module name (here, _interp) has to be the same as
                    # the basename (ie prefix) of the full filename.
                    '_covariance_matrix',
                    [ os.path.join('admixturebayes', '_covariance_matrix.pyx')],
                    include_dirs=[np.get_include()],
                    extra_compile_args = compile_args)
                ]
            )
else:
    extensions = [
        Extension(
            '_covariance_matrix', 
            [os.path.join('admixturebayes',  '_covariance_matrix.o')], 
        )
    ]

setup(
        name='AdmixtureBayes',
        version='0.3',
        license='Creative Commons Attribution-Noncommercial-Share Alike license',
        long_description=open('README.md').read(),
        ext_modules=extensions, 
        zip_safe=False, 
        
        install_requires=[
            "numpy>=1.11.0",
            "scipy>=0.17.0", 
            "matplotlib", 
            "pandas", 
            "pathos", 
	    "graphviz",
            "pygraphviz"
        ],
        
        packages=find_packages(),
        
        entry_points={'console_scripts': [
            'AdmixtureBayes = admixturebayes.__main__:main'
        ]},
        
        package_data={b'': ['admixturebayes/*.R']},
        include_package_data=True,
        )
        

