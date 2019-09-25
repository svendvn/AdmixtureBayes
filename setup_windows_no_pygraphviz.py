from setuptools import setup, find_packages

setup(
        name='AdmixtureBayes',
        version='0.2',
        license='Creative Commons Attribution-Noncommercial-Share Alike license',
        long_description=open('README.md').read(),
        zip_safe=False, 
        
        install_requires=[
            "numpy>=1.11.0",
            "scipy>=0.17.0", 
            "matplotlib", 
            "pandas", 
            "pathos", 
        ],
        
        packages=find_packages(),
        
        entry_points={'console_scripts': [
            'AdmixtureBayes = admixturebayes.__main__:main'
        ]},
        
        package_data={b'': ['admixturebayes/*.R']},
        include_package_data=True,
        )
        

