import os
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='BioLectorPy',
    version='0.2.1',    
    description='A Python package for plotting BioLector results',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/sulheim/biolector',
    author='Snorre sulheim',
    author_email='ssulheim@gmail.com',
    license='GPL-3.0',
    package_dir={"": str("src")},
    packages=find_packages(where="src", include=["BioLector"]),
    install_requires=['pandas',                     
                      'seaborn',
                      'matplotlib',
                      'pygam'],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',  
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',      
        'Programming Language :: Python :: 3',
    ],
)
