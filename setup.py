#!/usr/bin/env python
try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit("error: install setuptools")

setup(
    name='BioRanges',
    version=0.1,
    author='Vince Buffalo',
    author_email='vsbuffalo@gmail.com',
    packages=['BioRanges'],
    package_dir={'BioRanges': 'BioRanges'},
    url='http://github.com/vsbuffalo/BioRanges',
    license='GPL 2.0',
    description='A tiny lightweight ranges library, for BLAST primarily',
    requires=["BioPython"]
    )   
