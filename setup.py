import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='BioLectorPy',
    version='0.1.0',    
    description='A Python package for plotting BioLector results',
    long_description = read('README.md'),
    url='https://github.com/sulheim/biolector',
    author='Snorre sulheim',
    author_email='ssulheim@gmail.com',
    license='GPL-3.0',
    packages=['BioLector'],
    install_requires=['pandas',
                      'numpy',                     
                      'seaborn',
                      'matplotlib',
                      'pathlib'],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',  
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',      
        'Programming Language :: Python :: 3',
    ],
)
